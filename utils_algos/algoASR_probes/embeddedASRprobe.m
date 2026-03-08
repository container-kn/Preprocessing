classdef embeddedASRprobe < vanillaASRprobe
    % embeddedASRprobe - Diagnostic probe for embeddedASR.
    %
    %   Extends vanillaASRprobe (which itself extends asrProbe) with metrics
    %   specific to the time-delay embedding approach:
    %
    %   ── GROUP 1: Timing breakdown ────────────────────────────────────────
    %   Populated by calling finaliseTimings(asrObj) after all process() calls.
    %   Shows how much wall time each stage costs, including the two stages
    %   unique to embeddedASR that do not exist in vanillaASR:
    %
    %     t_embed        — trajectory matrix construction (embed_input stage)
    %     t_reconstruct  — anti-diagonal averaging back to 1D (reconstruct stage)
    %
    %   Shared ASR stages (also present in vanillaASR for direct comparison):
    %     t_spectralShaping, t_localCovariance, t_subspaceEig,
    %     t_riemannianDistance, t_reconstructionMatrix, t_calibration
    %
    %   Each timing field is a struct with fields: mean, std, min, max, n (seconds).
    %   t_calibration is a single scalar (one-shot).
    %
    %   ── GROUP 2: Embedding geometry (per process() call) ─────────────────
    %   Recorded once per process() call via recordEmbedding():
    %
    %     embedRatio       — L / N_orig: fraction of samples retained after
    %                        embedding (= 1 - (M-1)*tau/N). Below ~0.8 means
    %                        significant boundary loss.
    %     embeddedLength   — L: actual trajectory matrix column count
    %     reconstructCoverage  — fraction of output samples with counts > 0
    %                            (sanity check: should be 1.0 always)
    %     reconstructMeanCount — mean anti-diagonal overlap per output sample
    %                            (max = embDim when lag=1; less at boundaries)
    %     reconstructMaxCount  — maximum overlap at any single output sample
    %
    %   ── GROUP 3: Per-window subspace metrics (inherited) ─────────────────
    %   All asrProbe base metrics work in the embedded space:
    %     traceCov, condCov, spectralEntropy, riemannDrift,
    %     eigInflationMax/TopK, dirInflation*, energyRatio, normR, rankR,
    %     trivialFlag, logDetDrift
    %
    %   vanillaASRprobe additions (also inherited):
    %     subspaceAngleMax, subspaceAngleMean  — angle between V_t and V_0
    %     thresholdDrift                       — ||T_t - T_0|| / ||T_0||
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       obj.probeRaw   = embeddedASRprobe(embDim);
    %       obj.probeClean = embeddedASRprobe(embDim);
    %       obj.probeRaw.setBaseline(obj);       % after calibrate()
    %       obj.probeClean.setBaseline(obj);
    %
    %       % Inside process() — per ASR window (same call as vanillaASR):
    %       obj.probeRaw.update(Ct, Vw_desc, lambda_desc, threshold_dir_local, ...
    %                           R, trivial, riem_dist, obj.T_thresh);
    %
    %       % Inside process() — once per call:
    %       obj.probeRaw.recordEmbedding(N_orig, L_embedded, counts_vec);
    %
    %       % After all segments are processed:
    %       obj.probeRaw.finaliseTimings(obj);   % obj = embeddedASR instance
    %
    %       % Retrieve timing breakdown:
    %       obj.probeRaw.timingReport()
    %       obj.probeRaw.t_embed          % struct: mean/std/min/max/n
    %       obj.probeRaw.t_reconstruct    % struct: mean/std/min/max/n
    %
    %   ── setBaseline override ─────────────────────────────────────────────
    %   embeddedASR stores its mixing matrix as Mmix (not M) to avoid collision
    %   with the embedding dimension M. This probe's setBaseline() handles both
    %   naming conventions so it works with embeddedASR and vanillaASR alike.

    % ================================================================
    properties (Access = public)

        % ── Embedding geometry (per process() call) ──────────────────
        embedRatio           % [K x 1]  L / N_orig per call
        embeddedLength       % [K x 1]  L = trajectory columns per call
        reconstructCoverage  % [K x 1]  fraction of samples with count > 0
        reconstructMeanCount % [K x 1]  mean anti-diagonal overlap
        reconstructMaxCount  % [K x 1]  max anti-diagonal overlap

        % ── Timing structs (populated by finaliseTimings) ─────────────
        % Embedding-unique stages
        t_embed              % trajectory matrix construction time
        t_reconstruct        % anti-diagonal averaging time

        % Shared ASR stages
        t_spectralShaping    % IIR filter on embedded matrix
        t_localCovariance    % Ct estimation per window
        t_subspaceEig        % eig(Ct) per window
        t_riemannianDistance % AIRM per window
        t_reconstructionMatrix % pinv(…) per window
        t_calibration        % full calibrate() wall time (scalar struct)

    end % properties

    % ================================================================
    methods

        function obj = embeddedASRprobe(subspaceRank)
            % EMBEDDEDASRPROBE
            %   subspaceRank : number of virtual channels to track in subspace
            %                  metrics. Pass embDim (the embedding dimension).
            if nargin < 1, subspaceRank = []; end
            obj = obj@vanillaASRprobe(subspaceRank);
            obj.resetOwn();
        end

        % ── setBaseline ───────────────────────────────────────────────
        function setBaseline(obj, asrObj)
            % SETBASELINE - Capture calibration baseline from an embeddedASR
            %   (or vanillaASR) instance. Handles Mmix vs M naming.
            if ~obj.enabled, return; end

            % embeddedASR stores mixing matrix as Mmix to avoid name collision
            % with the embedding dimension. Build a temporary struct that
            % presents it as M so the parent setBaseline works unmodified.
            if isprop(asrObj, 'Mmix') && ~isempty(asrObj.Mmix)
                proxy   = struct();
                proxy.M = asrObj.Mmix;
                % T_thresh in embeddedASR corresponds to T in vanillaASR
                if isprop(asrObj, 'T_thresh') && ~isempty(asrObj.T_thresh)
                    proxy.T = asrObj.T_thresh;
                end
                setBaseline@vanillaASRprobe(obj, proxy);
            else
                % Fallback: vanillaASR-style object with .M property
                setBaseline@vanillaASRprobe(obj, asrObj);
            end
        end

        % ── recordEmbedding ───────────────────────────────────────────
        function recordEmbedding(obj, N_orig, L_embedded, counts_vec)
            % RECORDEMBEDDING - Capture embedding geometry for one process() call.
            %   N_orig      : length of the raw input signal [scalar]
            %   L_embedded  : number of trajectory matrix columns [scalar]
            %   counts_vec  : [1 x N_orig] overlap counts from reconstruct()
            %
            % Call once per process() call, after embed() and reconstruct().

            if ~obj.enabled, return; end

            obj.embeddedLength(end+1,1)  = L_embedded;
            obj.embedRatio(end+1,1)      = L_embedded / max(N_orig, 1);

            if nargin >= 4 && ~isempty(counts_vec)
                nz = counts_vec(counts_vec > 0);
                obj.reconstructCoverage(end+1,1)  = numel(nz) / max(N_orig, 1);
                obj.reconstructMeanCount(end+1,1) = mean(counts_vec);
                obj.reconstructMaxCount(end+1,1)  = max(counts_vec);
            else
                obj.reconstructCoverage(end+1,1)  = NaN;
                obj.reconstructMeanCount(end+1,1) = NaN;
                obj.reconstructMaxCount(end+1,1)  = NaN;
            end
        end

        % ── finaliseTimings ───────────────────────────────────────────
        function finaliseTimings(obj, asrObj)
            % FINALISETIMINGS - Pull timing summaries from the tProbe inside asrObj.
            %   Call once after all process() calls are complete.
            %   asrObj : the embeddedASR instance whose tProbe holds the records.

            if ~obj.enabled, return; end

            if isempty(asrObj) || ~isprop(asrObj, 'tProbe') || isempty(asrObj.tProbe)
                warning('embeddedASRprobe:noTProbe', ...
                    'asrObj.tProbe is empty — timing data unavailable.');
                return;
            end

            tp = asrObj.tProbe;

            obj.t_embed              = obj.safeStage(tp, 'embed_input');
            obj.t_reconstruct        = obj.safeStage(tp, 'reconstruct');
            obj.t_spectralShaping    = obj.safeStage(tp, 'spectral_shaping');
            obj.t_localCovariance    = obj.safeStage(tp, 'local_covariance');
            obj.t_subspaceEig        = obj.safeStage(tp, 'subspace_eig');
            obj.t_riemannianDistance = obj.safeStage(tp, 'riemannian_distance');
            obj.t_reconstructionMatrix = obj.safeStage(tp, 'reconstruction_matrix');
            obj.t_calibration        = obj.safeStage(tp, 'calibration');
        end

        % ── timingReport ──────────────────────────────────────────────
        function timingReport(obj)
            % TIMINGREPORT - Print a formatted timing breakdown to the console.
            %   Shows mean ± std and total time for each stage.
            %   Embedding-unique stages are highlighted with [E-ASR only].

            if ~obj.enabled
                fprintf('embeddedASRprobe: disabled.\n');
                return;
            end

            fprintf('\n══════════════════════════════════════════════════════\n');
            fprintf('  embeddedASRprobe — Timing Breakdown\n');
            fprintf('══════════════════════════════════════════════════════\n');

            stages = { ...
                't_calibration',         'Calibration (one-shot)',     false; ...
                't_embed',               'Trajectory matrix build',    true;  ...
                't_spectralShaping',     'Spectral shaping (IIR)',      false; ...
                't_localCovariance',     'Local covariance Ct',         false; ...
                't_subspaceEig',         'Eigendecomposition eig(Ct)', false; ...
                't_riemannianDistance',  'Riemannian distance (AIRM)',  false; ...
                't_reconstructionMatrix','Reconstruction matrix pinv',  false; ...
                't_reconstruct',         'Anti-diagonal averaging',     true;  ...
            };

            total_embed_time = 0;
            total_asr_time   = 0;

            for i = 1:size(stages, 1)
                fname  = stages{i,1};
                label  = stages{i,2};
                is_new = stages{i,3};

                s = obj.(fname);
                if isempty(s) || s.n == 0
                    fprintf('  %-32s  [no data]\n', label);
                    continue;
                end

                tag   = '';
                if is_new, tag = ' [E-ASR only]'; end

                total = s.mean * s.n;

                fprintf('  %-32s  mean=%7.4f ms  std=%7.4f ms  n=%5d  total=%7.3f s%s\n', ...
                    label, s.mean*1e3, s.std*1e3, s.n, total, tag);

                if is_new
                    total_embed_time = total_embed_time + total;
                else
                    total_asr_time = total_asr_time + total;
                end
            end

            fprintf('──────────────────────────────────────────────────────\n');
            fprintf('  Embedding overhead (E-ASR only) :  %.3f s\n', total_embed_time);
            fprintf('  ASR pipeline total               :  %.3f s\n', total_asr_time);
            if total_asr_time > 0
                fprintf('  Embedding / ASR ratio            :  %.1f%%\n', ...
                    100 * total_embed_time / (total_embed_time + total_asr_time));
            end
            fprintf('══════════════════════════════════════════════════════\n\n');
        end

        % ── embeddingReport ───────────────────────────────────────────
        function embeddingReport(obj)
            % EMBEDDINGREPORT - Print embedding geometry summary across all calls.

            if isempty(obj.embedRatio)
                fprintf('embeddedASRprobe: no embedding records yet.\n');
                return;
            end

            fprintf('\n══════════════════════════════════════════════════════\n');
            fprintf('  embeddedASRprobe — Embedding Geometry\n');
            fprintf('══════════════════════════════════════════════════════\n');
            fprintf('  process() calls recorded : %d\n',   numel(obj.embedRatio));
            fprintf('  Embed ratio  mean / min  : %.3f / %.3f\n', ...
                mean(obj.embedRatio), min(obj.embedRatio));
            fprintf('  Trajectory L mean / min  : %.0f / %.0f samples\n', ...
                mean(obj.embeddedLength), min(obj.embeddedLength));
            fprintf('  Reconstruct coverage     : %.4f (should be 1.0)\n', ...
                mean(obj.reconstructCoverage, 'omitnan'));
            fprintf('  Overlap mean / max       : %.2f / %.2f\n', ...
                mean(obj.reconstructMeanCount, 'omitnan'), ...
                mean(obj.reconstructMaxCount, 'omitnan'));
            fprintf('══════════════════════════════════════════════════════\n\n');
        end

        % ── reset ─────────────────────────────────────────────────────
        function reset(obj)
            reset@vanillaASRprobe(obj);
            obj.resetOwn();
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function resetOwn(obj)
            obj.embedRatio           = zeros(0,1);
            obj.embeddedLength       = zeros(0,1);
            obj.reconstructCoverage  = zeros(0,1);
            obj.reconstructMeanCount = zeros(0,1);
            obj.reconstructMaxCount  = zeros(0,1);
            obj.t_embed              = [];
            obj.t_reconstruct        = [];
            obj.t_spectralShaping    = [];
            obj.t_localCovariance    = [];
            obj.t_subspaceEig        = [];
            obj.t_riemannianDistance = [];
            obj.t_reconstructionMatrix = [];
            obj.t_calibration        = [];
        end

        function s = safeStage(~, tp, stageName)
            % SAFESTAGE - Summarize a tProbe stage; return empty struct on miss.
            try
                s = tp.summarize(stageName);
            catch
                s = struct('mean', NaN, 'std', NaN, 'min', NaN, 'max', NaN, 'n', 0);
            end
        end

    end % private methods

end
