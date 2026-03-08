classdef embeddedASR < handle
    % embeddedASR - Artifact Subspace Reconstruction via time-delay embedding.
    %
    %   Applies standard ASR in a delay-coordinate (trajectory matrix) space.
    %   A single-channel EEG signal is lifted to M virtual channels using
    %   Takens-style embedding:
    %
    %       E[m, l] = signal[l + (m-1)*tau]    m = 1..M,  l = 1..L
    %
    %   This M x L matrix is treated as M-channel EEG and processed by the
    %   vanillaASR subspace-reconstruction pipeline.  The cleaned matrix is
    %   inverted back to 1D via anti-diagonal (Hankel) averaging.
    %
    %   The key advantage: artefacts that are low-rank in delay-coordinate
    %   space (e.g. blinks, which are stereotyped waveforms) are separated
    %   from neural background even when the raw signal is single-channel.
    %
    %   ── Embedding parameters ────────────────────────────────────────────
    %     embDim    M  — number of delay coordinates.
    %                    Default: ceil(srate / fL) where fL = 1 Hz.
    %                    Rule of thumb: at least one full period of the slowest
    %                    frequency of interest.
    %     lag       τ  — time lag between successive coordinates (samples).
    %                    Default: 1 (consecutive samples).
    %
    %   ── Input / output convention ────────────────────────────────────────
    %     calibrate(X) : X is [1 x S] or [S x 1] single-channel EEG.
    %     process(X)   : same — returns [1 x S] cleaned single-channel EEG.
    %     modifiedMask : [1 x S] logical — true where reconstruction fired.
    %
    %   ── Constructor signature (matches Experimenting.m harness) ─────────
    %       obj = embeddedASR(srate, params)
    %       obj = embeddedASR(srate, params, tProbe)
    %
    %   ── Experimenting.m case ─────────────────────────────────────────────
    %       case {'e-asr', 'e_asr', 'embeddedasr'}
    %           obj.asr = embeddedASR(fs, params, tp);
    %
    %   ── Parameters (set via params struct or directly on object) ────────
    %     embDim             M  embedding dimension  (default: ceil(srate/1))
    %     lag                τ  time lag in samples  (default: 1)
    %     cutoff             k  SD threshold         (default: 15)
    %     blocksize             covariance block size (default: 10)
    %     window_len            processing window (s) (default: 0.5)
    %     window_overlap        window overlap frac   (default: 0.66)
    %     maxdims               max dims to reject    (default: 0.66)
    %     lookahead             lookahead buffer (s)  (default: 0.25)
    %     stepsize              hop between updates   (default: 32)
    %
    %   ── Reference ────────────────────────────────────────────────────────
    %     Based on embeddedASR_working.m (Vishnu K., IIT Guwahati).
    %     Embedding follows Takens (1981); ASR follows Mullen et al. (2015).

    % ================================================================
    properties (Access = public)

        % ── Core ─────────────────────────────────────────────────────
        srate
        nchans      = 1      % physical channels (always 1 for single-channel input)
        nsamples    = 0
        modifiedMask

        % ── Embedding parameters ─────────────────────────────────────
        embDim               % M: number of delay coordinates
        lag         = 1      % τ: time lag between coordinates (samples)

        % ── ASR parameters ───────────────────────────────────────────
        cutoff              = 15     % SD threshold k
        blocksize           = 10     % covariance block size
        window_len          = 0.5    % processing window length (s)
        window_overlap      = 0.66   % overlap fraction
        max_dropout_fraction = 0.1
        min_clean_fraction   = 0.25
        maxdims             = 0.66   % max fraction of virtual channels to reject
        lookahead           = 0.25   % lookahead buffer (s)
        stepsize            = 32     % samples between successive updates
        epsi                = 1e-12

        % ── Calibration state ────────────────────────────────────────
        % NOTE: Mmix is the mixing matrix (sqrt of embedded covariance).
        %       Named Mmix to avoid shadowing the embedding dimension M.
        Mmix                 % sqrt of baseline covariance  [embDim x embDim]
        T_thresh             % threshold operator matrix    [embDim x embDim]
        V_eig                % eigenvectors of Mmix         [embDim x embDim]
        threshold_dir        % per-component variance thresholds
        C0_inv_sqrt          % C0^{-1/2} for AIRM computation

        % ── Spectral filter ──────────────────────────────────────────
        B                    % Yule-Walker numerator
        A                    % Yule-Walker denominator

        % ── Runtime buffers ──────────────────────────────────────────
        carry       = []     % lookahead boundary buffer [embDim x P]
        iir_state   = []     % IIR filter state
        last_R      = []     % reconstruction matrix from previous step
        last_trivial = true  % trivial flag from previous step
        blockCount  = 0

        % ── Telemetry ────────────────────────────────────────────────
        tProbe               % timeProbe instance
        probeRaw             % vanillaASRprobe — diagnostics on embedded signal
        probeClean           % vanillaASRprobe — diagnostics on cleaned embedded signal

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = embeddedASR(srate, params, tProbe)
            % EMBEDDEDASR
            %   srate  : sampling rate Hz (required)
            %   params : struct of parameter overrides (optional)
            %   tProbe : timeProbe instance (optional — created if absent)

            if nargin < 1 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('embeddedASR:invalidSrate', 'srate must be a positive numeric scalar.');
            end
            obj.srate = srate;

            % tProbe
            if nargin >= 3 && ~isempty(tProbe)
                obj.tProbe = tProbe;
            else
                obj.tProbe = timeProbe();
            end

            % Default embedding dimension: one full period of lowest freq (1 Hz)
            obj.embDim = ceil(srate / 1);

            % Apply parameter overrides
            if nargin >= 2 && isstruct(params)
                f = fieldnames(params);
                for i = 1:numel(f)
                    if isprop(obj, f{i})
                        obj.(f{i}) = params.(f{i});
                    end
                end
                % Legacy field name used by configExperiment: embeddingDimension
                if isfield(params, 'embeddingDimension') && ~isempty(params.embeddingDimension)
                    obj.embDim = params.embeddingDimension;
                end
            end

            obj.embDim = max(2, round(obj.embDim));   % must be at least 2

            [obj.B, obj.A] = embeddedASR.buildSpectralFilter(srate);
            obj.modifiedMask = false(1, 0);
        end

        % ── calibrate ────────────────────────────────────────────────
        function calibrate(obj, X)
            % CALIBRATE - Build clean model from single-channel baseline EEG.
            %   X : [1 x S] or [S x 1] single-channel calibration signal.

            if isempty(X) || ~isnumeric(X)
                error('embeddedASR:calibrate', 'Calibration data must be non-empty numeric.');
            end
            X = embeddedASR.toRowVec(X);

            fprintf('--- [embeddedASR] Calibrating (embDim=%d, lag=%d) ---\n', ...
                obj.embDim, obj.lag);
            obj.tProbe.start('calibration');

            % 1. Embed calibration signal → trajectory matrix
            E_calib = obj.embed(X);   % [embDim x L_calib]

            % 2. Spectral shaping on embedded matrix
            obj.tProbe.start('prefilter');
            E_filt = filter(obj.B, obj.A, double(E_calib), [], 2);
            obj.tProbe.stop('prefilter');

            % 3. Block-averaged covariance of embedded calibration signal
            M = obj.embDim;
            U = embeddedASR.blockCovariance(E_filt, obj.blocksize);
            U = (U + U') / 2;
            U = U + obj.epsi * eye(M);

            % 4. Mixing matrix (sqrt of baseline covariance)
            obj.tProbe.start('subspaceM');
            obj.Mmix = sqrtm(real(U));
            obj.tProbe.stop('subspaceM');

            % 5. C0^{-1/2} for AIRM
            [V0, D0]  = eig(obj.Mmix * obj.Mmix');
            lambda0   = diag(D0);
            lambda0(lambda0 < obj.epsi) = obj.epsi;
            obj.C0_inv_sqrt = V0 * diag(1 ./ sqrt(lambda0)) * V0';

            % 6. Eigenvectors for threshold projection
            obj.tProbe.start('eigCal');
            [obj.V_eig, ~] = eig(obj.Mmix);
            obj.tProbe.stop('eigCal');

            % 7. Threshold estimation
            [obj.T_thresh, obj.threshold_dir] = obj.getThresholdVector(E_filt);

            % 8. Probes — embeddedASRprobe tracks embedding geometry + timing
            obj.probeRaw   = embeddedASRprobe(obj.embDim);
            obj.probeClean = embeddedASRprobe(obj.embDim);
            obj.probeRaw.setBaseline(obj);
            obj.probeClean.setBaseline(obj);

            % 9. Reset runtime state
            obj.carry        = [];
            obj.last_R       = [];
            obj.last_trivial = true;
            obj.iir_state    = [];
            obj.blockCount   = 0;

            obj.tProbe.stop('calibration');
            fprintf('--- [embeddedASR] Calibration complete ---\n');
        end

        % ── process ──────────────────────────────────────────────────
        function out = process(obj, data)
            % PROCESS - Apply embedded ASR to a single-channel EEG segment.
            %   data : [1 x S] or [S x 1] single-channel EEG.
            %   out  : [1 x S] cleaned single-channel EEG.

            if isempty(obj.Mmix)
                error('embeddedASR:notCalibrated', 'Run calibrate() before process().');
            end
            if isempty(data) || ~isnumeric(data)
                error('embeddedASR:invalidInput', 'Data must be non-empty numeric.');
            end

            data     = embeddedASR.toRowVec(data);
            N_orig   = length(data);
            M        = obj.embDim;
            C        = M;    % number of "virtual" channels
            P        = round(obj.lookahead * obj.srate);

            % ── 1. Embed input signal ────────────────────────────────
            % Prepend carry (or mirror-reflect on first call) before embedding
            % so that the boundary is handled in the ORIGINAL signal domain.
            if isempty(obj.carry)
                % Mirror-reflect first P samples for causal boundary
                mirror_idx = 1 + mod(((P+1):-1:2) - 1, N_orig);
                pre        = 2 * data(1) - data(mirror_idx);   % [1 x P]
                data_pad   = [pre, double(data)];
            else
                % carry is [1 x P] from the previous call
                data_pad = [obj.carry, double(data)];
            end
            data_pad(~isfinite(data_pad)) = 0;

            % Extend modifiedMask for incoming samples
            obj.modifiedMask(obj.nsamples + (1:N_orig)) = false;

            % Embed the padded signal — timed separately (embedding-unique cost)
            obj.tProbe.start('embed_input');
            E = obj.embed(data_pad);   % [M x L_pad]
            obj.tProbe.stop('embed_input');

            % ── 2. Spectral shaping ──────────────────────────────────
            obj.tProbe.start('spectral_shaping');
            [X, obj.iir_state] = filter(obj.B, obj.A, double(E), obj.iir_state, 2);
            obj.tProbe.stop('spectral_shaping');

            % ── 3. ASR processing window loop ────────────────────────
            % Identical to vanillaASR.process() but operating on embedded space.
            N_win    = round(obj.window_len * obj.srate);
            update_at = obj.stepsize : obj.stepsize : size(X, 2);

            if isempty(obj.last_R)
                obj.last_R = eye(C);
                update_at  = [1, update_at];
            end

            last_n = 0;

            for j = 1:length(update_at)
                obj.blockCount = obj.blockCount + 1;
                n   = update_at(j);
                idx = max(1, n - N_win + 1) : n;
                Xw  = X(:, idx);

                % Local covariance
                obj.tProbe.start('local_covariance');
                Ct = (Xw * Xw') / size(Xw, 2);
                Ct = (Ct + Ct') / 2 + 1e-10 * eye(C);
                obj.tProbe.stop('local_covariance');

                % Eigendecomposition
                obj.tProbe.start('subspace_eig');
                [Vw, Dmat]          = eig(Ct);
                lambda_raw          = diag(Dmat);
                lambda_raw(lambda_raw <= 0) = eps;
                [lambda_asc, asc_ord]  = sort(lambda_raw, 'ascend');
                Vw_asc              = Vw(:, asc_ord);
                [lambda_desc, desc_ord] = sort(lambda_raw, 'descend');
                Vw_desc             = Vw(:, desc_ord);
                obj.tProbe.stop('subspace_eig');

                % AIRM Riemannian distance
                obj.tProbe.start('riemannian_distance');
                delta     = obj.C0_inv_sqrt * Ct * obj.C0_inv_sqrt;
                delta     = (delta + delta') / 2;
                lam_delta = eig(delta);
                riem_dist = sqrt(sum(log(max(lam_delta, 1e-12)).^2));
                obj.tProbe.stop('riemannian_distance');

                % Threshold evaluation
                threshold_dir_local = sum((obj.T_thresh * Vw_asc).^2);
                keep    = (lambda_asc < threshold_dir_local') | ...
                          ((1:C)' <= C - round(obj.maxdims * C));
                trivial = all(keep);

                % Reconstruction matrix
                obj.tProbe.start('reconstruction_matrix');
                if ~trivial
                    R = real(obj.Mmix * pinv(bsxfun(@times, keep', Vw_asc' * obj.Mmix)) * Vw_asc');
                else
                    R = eye(C);
                end
                obj.tProbe.stop('reconstruction_matrix', ...
                    struct('trivial', trivial, 'dims_rejected', sum(~keep)));

                % Signal blending + probe updates
                if ~trivial || ~obj.last_trivial
                    subrange    = (last_n + 1) : n;
                    raw_segment = E(:, subrange);   % blend in EMBEDDED domain

                    % Map subrange in embedded domain back to original-signal samples.
                    % The embedded column l corresponds to original sample l (1-indexed
                    % within data_pad). Strip the P-sample prefix to get mask indices
                    % in the output frame.
                    mask_idx = obj.nsamples + subrange - P;
                    mask_idx = mask_idx(mask_idx >= 1 & mask_idx <= obj.nsamples + N_orig);
                    obj.modifiedMask(mask_idx) = true;

                    blend = (1 - cos(pi * (1:length(subrange)) / length(subrange))) / 2;
                    cleaned_seg = bsxfun(@times, blend,   R          * raw_segment) + ...
                                  bsxfun(@times, 1-blend, obj.last_R * raw_segment);
                    E(:, subrange) = cleaned_seg;

                    % Probe updates on embedded signal
                    obj.probeRaw.update(Ct, Vw_desc, lambda_desc, threshold_dir_local, ...
                                        R, trivial, riem_dist, obj.T_thresh);

                    Ct_c = (cleaned_seg * cleaned_seg') / size(cleaned_seg, 2);
                    Ct_c = (Ct_c + Ct_c') / 2 + 1e-10 * eye(C);
                    [Vc, Dc] = eig(Ct_c);
                    lam_c    = diag(Dc);
                    [lam_cd, cord] = sort(lam_c, 'descend');
                    obj.probeClean.update(Ct_c, Vc(:,cord), lam_cd, threshold_dir_local, ...
                                          R, trivial, riem_dist, obj.T_thresh);
                end

                obj.last_R       = R;
                obj.last_trivial = trivial;
                last_n           = n;
            end

            % ── 4. Update carry buffer (in original signal domain) ───
            % Store last P samples of the unembedded padded signal for next call.
            obj.carry = data_pad(end-P+1 : end);   % [1 x P]

            % ── 5. Reconstruct 1D from cleaned embedded matrix ───────
            % E is [M x L_pad] in the padded frame. We need the columns that
            % correspond to original (non-carry) output samples P+1 : N_orig+P.
            % Then anti-diagonal average maps [M x L_out] back to [1 x N_orig].
            L_pad = size(E, 2);

            % Columns of E that start within the output window:
            % Column l maps to original-padded sample l (lag=1 case) or l to
            % l+(M-1)*tau. To get output sample range [P+1, N_orig+P] we need
            % the columns whose FIRST row (l itself) falls in that range.
            col_start = P + 1;
            col_end   = min(L_pad, N_orig + P);
            E_out     = E(:, col_start : col_end);   % [M x N_orig] approximately

            % Reconstruct — timed separately (embedding-unique cost)
            obj.tProbe.start('reconstruct');
            [out, counts] = obj.reconstruct(E_out, N_orig);   % [1 x N_orig]
            obj.tProbe.stop('reconstruct');

            % Record embedding geometry in probes (once per process() call)
            L_embedded = size(E, 2);
            obj.probeRaw.recordEmbedding(N_orig, L_embedded, counts);
            obj.probeClean.recordEmbedding(N_orig, L_embedded, counts);

            obj.nsamples = obj.nsamples + N_orig;
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        % ── embed ─────────────────────────────────────────────────────
        function E = embed(obj, signal)
            % EMBED - Build trajectory matrix from 1D signal.
            %   signal : [1 x N] row vector
            %   E      : [embDim x L] where L = N - (embDim-1)*lag
            N   = length(signal);
            M   = obj.embDim;
            tau = obj.lag;
            L   = N - (M - 1) * tau;
            if L <= 0
                error('embeddedASR:tooShort', ...
                    'Signal length (%d) too short for embDim=%d, lag=%d (need >= %d).', ...
                    N, M, tau, (M-1)*tau + 1);
            end
            % Vectorized index matrix: idx(m, l) = l + (m-1)*tau
            idx = (1:L) + (0:M-1)' * tau;   % [M x L]
            E   = signal(idx);               % [M x L]
        end

        % ── reconstruct ───────────────────────────────────────────────
        function [out, counts] = reconstruct(obj, E, N_out)
            % RECONSTRUCT - Anti-diagonal (Hankel) averaging back to 1D.
            %   E     : [M x L] cleaned trajectory matrix
            %   N_out : desired output length
            %   out   : [1 x N_out]
            [M, L] = size(E);
            tau    = obj.lag;

            out    = zeros(1, N_out);
            counts = zeros(1, N_out);

            for m = 1:M
                idx_samples = (1:L) + (m - 1) * tau;
                % Only accumulate where we have valid output indices
                valid = idx_samples >= 1 & idx_samples <= N_out;
                out(idx_samples(valid))    = out(idx_samples(valid))    + E(m, valid);
                counts(idx_samples(valid)) = counts(idx_samples(valid)) + 1;
            end

            nz        = counts > 0;
            out(nz)   = out(nz) ./ counts(nz);
        end

        % ── getThresholdVector ────────────────────────────────────────
        function [T_out, threshold_dir] = getThresholdVector(obj, X)
            % Estimate per-component amplitude thresholds from calibration.
            % Mirrors vanillaASR.getThresholdVector exactly.
            [C, S] = size(X);
            N      = round(obj.window_len * obj.srate);
            Xp     = abs(obj.V_eig' * X);

            rms_vals = zeros(C, max(1, floor(S / N)));
            for k = 1:size(rms_vals, 2)
                rms_vals(:, k) = sqrt(mean(Xp(:, (k-1)*N + (1:N)).^2, 2));
            end

            T_vec = zeros(C, 1);
            for c = 1:C
                v = rms_vals(c, isfinite(rms_vals(c,:)));
                if isempty(v), T_vec(c) = obj.cutoff; continue; end
                try
                    [mu_c, sig_c] = fit_eeg_distribution(v, ...
                        obj.min_clean_fraction, obj.max_dropout_fraction);
                    T_vec(c) = mu_c + obj.cutoff * sig_c;
                catch
                    T_vec(c) = median(v) + obj.cutoff * mad(v, 1) * 1.4826;
                end
            end

            T_out         = diag(T_vec) * obj.V_eig';
            threshold_dir = sum((T_out * obj.V_eig).^2);
        end

    end % private methods

    % ================================================================
    methods (Static)

        function [B, A] = buildSpectralFilter(srate)
            % Yule-Walker IIR spectral shaping filter — identical to vanillaASR.
            nyq      = srate / 2;
            freqvals = [0 2 3 13 16 40 min(80, nyq-1)] * 2 / srate;
            freqvals = [freqvals, 1];
            amps     = [3 0.75 0.33 0.33 1 1 3 3];
            [B, A]   = yulewalk(8, freqvals, amps);
        end

        function U = blockCovariance(X, blocksize)
            % Block-averaged covariance — identical to vanillaASR.blockCovariance.
            [~, S] = size(X);
            U      = zeros(size(X, 1));
            n_used = 0;
            for k = 1:blocksize
                idx = k : blocksize : S;
                if isempty(idx), continue; end
                Xk     = X(:, idx);
                U      = U + (Xk * Xk');
                n_used = n_used + size(Xk, 2);
            end
            if n_used == 0
                error('embeddedASR:blockCovariance', 'No samples for covariance estimation.');
            end
            U = U / n_used;
        end

        function x = toRowVec(x)
            % Ensure input is a row vector [1 x N].
            x = x(:)';
        end

    end % static methods

end
