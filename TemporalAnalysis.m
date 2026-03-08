classdef TemporalAnalysis < handle
    % TemporalAnalysis - Windowed signal and subspace statistics over time.
    %
    %   Answers the question SignalStatistics cannot: WHERE in the recording
    %   did ASR introduce distortion, and how did signal geometry change
    %   window by window?
    %
    %   Computes two groups of statistics in a sliding window over each segment:
    %
    %   ── Signal group ─────────────────────────────────────────────────────
    %   Per window, mean across selected channels:
    %     rrmse        rms(clean-raw) / rms(raw)        reconstruction error
    %     corr         Pearson corr(raw, clean)          signal preservation
    %     snrDB        10*log10(power_raw / power_diff)  SNR per channel, then mean
    %     rms_raw/clean  RMS amplitude before/after
    %     std_raw/clean  standard deviation before/after
    %     mean_raw/clean mean amplitude before/after
    %
    %   ── Subspace group ───────────────────────────────────────────────────
    %   Per window, using top-k eigenvectors of per-window covariance matrices:
    %     angle_mean     mean principal angle (deg) between subspaces of Cx, Cy
    %     angle_max      max  principal angle (deg)
    %     angle_all      {1xW} full k-vector per window (stored for distribution)
    %     energy_ratio   trace(Cy) / trace(Cx)          — 1 = no energy change
    %     energy_removed 1 - energy_ratio                — fraction removed
    %     cov_fro_norm   ||Cy - Cx||_F / ||Cx||_F       — covariance distortion
    %     eigcorr        corr(eigenvalue spectrum Cx, Cy) — spectral similarity
    %
    %   ── Energy compaction group ──────────────────────────────────────────
    %   Per window, how compactly is energy organised in raw vs clean?
    %   (from compute_energy_compaction_time.m)
    %     kCompact_raw    dims needed for compactionThresh% energy in raw
    %     kCompact_clean  dims needed for compactionThresh% energy in cleaned
    %     dimReduction    kCompact_raw - kCompact_clean
    %                     >0: ASR compacted energy into fewer dims (good if
    %                         artifacts were high-rank)
    %                     <0: ASR spread energy (potentially over-correcting)
    %     energyCurve_raw   {1xW} full cumulative energy curve, raw  [C x 1]
    %     energyCurve_clean {1xW} same, cleaned
    %
    %   ── Calibration drift group ──────────────────────────────────────────
    %   Per window, how far is the current window from the calibration baseline?
    %   Reference subspace V_ref computed from raw.calibration once.
    %   (from compute_energy_capture_time.m + compute_subspace_angle_time.m)
    %
    %     energyCapture_raw    trace(V_ref' * Cx * V_ref) / trace(Cx)
    %                          fraction of raw window energy in calibration subspace
    %     energyCapture_clean  same for cleaned
    %     energyCaptureGain    energyCapture_clean - energyCapture_raw
    %                          >0: ASR pushed signal TOWARD calibration subspace
    %     angleFromRef_raw     mean principal angle between raw window and V_ref (deg)
    %     angleFromRef_clean   same for cleaned
    %     angleFromRefGain     angleFromRef_raw - angleFromRef_clean
    %                          >0: ASR moved signal CLOSER to calibration baseline
    %
    %   ── Energy compaction group ──────────────────────────────────────────
    %   Per window, how compactly is energy organised in raw vs clean?
    %   (implements compute_energy_compaction_time.m)
    %
    %     kCompact_raw    number of dims needed for compactionThresh% energy, raw
    %     kCompact_clean  same for cleaned signal
    %     dimReduction    kCompact_raw - kCompact_clean
    %                     >0: ASR compacted energy into fewer dims (desired — artifacts
    %                         are typically high-rank, brain activity is low-rank)
    %                     <0: ASR spread energy (over-correction signal)
    %     energyCurve_raw   {1xW} cumulative energy fractions [C x 1] each, raw
    %     energyCurve_clean {1xW} same for cleaned
    %
    %   ── Calibration drift group ──────────────────────────────────────────
    %   Per window vs a fixed reference subspace V_ref derived from the
    %   calibration segment. Answers: does cleaned signal stay within
    %   the baseline subspace better than the raw signal does?
    %   (implements compute_energy_capture_time.m + compute_subspace_angle_time.m
    %    + compute_state_covariance_rms.m for V_ref estimation)
    %
    %     energyCapture_raw    trace(V_ref' * Cx * V_ref) / trace(Cx)
    %                          fraction of raw window energy in calibration subspace
    %     energyCapture_clean  same for cleaned
    %     energyCaptureGain    energyCapture_clean - energyCapture_raw
    %                          >0: ASR pushed signal toward calibration (correct)
    %                          <0: ASR removed calibration-subspace energy (bad)
    %     angleFromRef_raw     mean principal angle (deg): raw window vs V_ref
    %     angleFromRef_clean   same for cleaned window
    %     angleFromRefGain     angleFromRef_raw - angleFromRef_clean
    %                          >0: ASR moved signal closer to calibration (correct)
    %
    %   V_ref is computed once from raw.calibration before any combo loop.
    %   If useRobustCov=true (default) and block_geometric_median (EEGLAB) is
    %   available, uses the ASR-compatible robust block geometric median
    %   covariance estimator. Falls back to simple cov() automatically.
    %
    %   Time-averaged scalars exposed as summary fields for paramTable/getSummary.
    %
    %   ── Results structure ────────────────────────────────────────────────
    %   results(c)
    %     .parameters
    %     .segment.(closed|open)
    %         .time             [1 x W]  window centre times (sec)
    %         .signal           struct of [1 x W] arrays
    %         .subspace         struct of [1 x W] arrays
    %         .energyCompaction struct of [1 x W] arrays + {1 x W} cell curves
    %         .calibrationDrift struct of [1 x W] arrays
    %         .summary
    %             .signal.(field)_{mean|std|max|p95}
    %             .subspace.(field)_{mean|std|max|p95}
    %             .energyCompaction.(field)_{mean|std|max|p95}
    %             .calibrationDrift.(field)_{mean|std|max|p95}
    %
    %   ── getSummary interface ─────────────────────────────────────────────
    %   fieldName = 'group.field_aggregation'  e.g.:
    %     'signal.rrmse_mean'                     mean RRMSE over windows
    %     'signal.corr_mean'                      mean correlation over windows
    %     'signal.rrmse_p95'                      95th percentile RRMSE
    %     'subspace.angle_mean_mean'              mean principal angle (raw vs clean)
    %     'subspace.cov_fro_norm_mean'            mean covariance distortion
    %     'subspace.energy_removed_max'           peak energy removal
    %     'energyCompaction.kCompact_raw_mean'    mean dims for thresh% energy, raw
    %     'energyCompaction.kCompact_clean_mean'  mean dims for thresh% energy, clean
    %     'energyCompaction.dimReduction_mean'    mean dimensionality reduction
    %     'calibrationDrift.energyCaptureGain_mean'  mean energy capture improvement
    %     'calibrationDrift.angleFromRefGain_mean'   mean angle-from-ref improvement
    %     'subspace.cov_fro_norm_mean' mean covariance distortion
    %     'subspace.energy_removed_max' peak energy removal
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       ta = TemporalAnalysis(raw, subject_results, fs, subjectID, algo);
    %       ta.compute('windowSec', 2, 'hopSec', 1, 'subspaceRank', 10);
    %       ta.getSummary('signal.rrmse_mean',        'closed')
    %       ta.getSummary('subspace.angle_mean_mean', 'closed')
    %       ta.plotTimeSeries(1, 'signal.rrmse',      'closed')
    %       ta.plotComparison('signal.corr',          'closed')
    %       T = ta.paramTable();
    %
    %   ── Relationship to other classes ────────────────────────────────────
    %   SignalStatistics   — scalar per segment; TemporalAnalysis = time-resolved
    %   QuietRegionAnalysis — metrics on quiet samples; this is ALL samples, windowed
    %   probeAnalysis      — subspaceAngle from ASR internals; this is external,
    %                        computed independently from raw/clean signals
    %   SpectralAnalysis   — frequency domain; this is time domain geometry

    % ================================================================
    properties (Access = public)

        raw               % struct: .closed/.open [C x T]
        subject_results   % 1xN struct array (.cleanClosed/.cleanOpen)
        fs
        subjectID
        algorithmName
        nCombos

        % Window config — may be set before compute()
        windowSec    = 2      % window length (sec)
        hopSec       = 1      % hop size (sec); [] = non-overlapping (= windowSec)
        subspaceRank = []     % top-k eigenvectors for subspace group; [] = min(10,nCh)
        segments     = {'closed', 'open'}

        % Energy compaction parameters
        % (compute_energy_compaction_time.m)
        compactionThresh = 0.9   % cumulative energy threshold for kCompact (0–1)

        % Calibration drift parameters
        % (compute_energy_capture_time.m + compute_subspace_angle_time.m)
        kRef          = []     % reference subspace rank; [] = same as kRank
        useRobustCov  = true   % use block_geometric_median (EEGLAB) if available
        blockSize     = 10     % blocksize for robust cov (compute_state_covariance_rms)

        % Optional streaming loader — set by ExperimentAnalysis when using
        % split-file layout.  Signature: loaderFcn(k) loads combo k into
        % obj.subject_results(k) and unloaderFcn(k) frees it afterward.
        loaderFcn   = []   % @(k) ana.loadCombo(k)
        unloaderFcn = []   % @(k) ana.unloadCombo(k)

        % Results — 1xN struct array (populated by compute())
        results

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = TemporalAnalysis(raw, subject_results, fs, subjectID, algorithmName)
            % TEMPORALANALYSIS
            %   raw            : struct with .closed/.open [C x T]
            %   subject_results: 1xN struct array (.cleanClosed/.cleanOpen)
            %   fs             : sampling rate Hz
            %   subjectID      : scalar  (default 0)
            %   algorithmName  : string  (default 'unknown')

            if nargin < 4, subjectID     = 0;         end
            if nargin < 5, algorithmName = 'unknown'; end

            obj.raw             = raw;
            obj.subject_results = subject_results;
            obj.fs              = fs;
            obj.subjectID       = subjectID;
            obj.algorithmName   = algorithmName;
            obj.nCombos         = numel(subject_results);
            obj.results         = struct();
        end

        % ── compute ──────────────────────────────────────────────────
        function compute(obj, varargin)
            % COMPUTE - Run windowed analysis across all combos and segments.
            %
            %   Name-value options:
            %     'windowSec'    2       Window length in seconds
            %     'hopSec'       1       Hop size in seconds ([] = non-overlapping)
            %     'subspaceRank' []      Top-k eigenvectors ([] = min(10,nCh))
            %     'channels'     []      Channel subset  ([] = all)

            p = inputParser;
            addParameter(p, 'windowSec',         obj.windowSec,         @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'hopSec',             obj.hopSec,            @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
            addParameter(p, 'overlap',            [],                    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0 && x < 1));
            addParameter(p, 'subspaceRank',       obj.subspaceRank,      @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
            addParameter(p, 'channels',           [],                    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
            addParameter(p, 'compactionThresh',   obj.compactionThresh,  @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
            addParameter(p, 'kRef',               obj.kRef,              @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
            addParameter(p, 'useRobustCov',       obj.useRobustCov,      @(x) islogical(x) || isnumeric(x));
            addParameter(p, 'blockSize',          obj.blockSize,         @(x) isnumeric(x) && isscalar(x) && x >= 1);
            parse(p, varargin{:});

            obj.windowSec       = p.Results.windowSec;
            obj.hopSec          = p.Results.hopSec;
            obj.subspaceRank    = p.Results.subspaceRank;
            obj.compactionThresh = p.Results.compactionThresh;
            obj.kRef            = p.Results.kRef;
            obj.useRobustCov    = logical(p.Results.useRobustCov);
            obj.blockSize       = p.Results.blockSize;
            chParam             = p.Results.channels;

            % overlap overrides hopSec if provided
            ov = p.Results.overlap;
            if ~isempty(ov)
                obj.hopSec = obj.windowSec * (1 - ov);
            end

            if isempty(obj.hopSec)
                obj.hopSec = obj.windowSec;   % non-overlapping
            end

            winSamp = round(obj.windowSec * obj.fs);
            hopSamp = round(obj.hopSec    * obj.fs);

            fprintf('TemporalAnalysis: S%d %s (%d combo(s), win=%.1fs hop=%.1fs)...\n', ...
                obj.subjectID, obj.algorithmName, obj.nCombos, ...
                obj.windowSec, obj.hopSec);

            % Raw signal map — same across all combos
            rawMap = struct('closed', double(obj.raw.closed), ...
                            'open',   double(obj.raw.open));

            nCh = size(obj.raw.closed, 1);
            if isempty(chParam)
                chList = 1:nCh;
            else
                chList = chParam(:)';
                chList = chList(chList >= 1 & chList <= nCh);
            end
            nChSel = numel(chList);

            % Subspace rank: min(k, nChSel) — enforce at compute time
            if isempty(obj.subspaceRank)
                kRank = min(10, nChSel);
            else
                kRank = min(obj.subspaceRank, nChSel);
            end

            % Reference subspace rank for calibration drift metrics
            if isempty(obj.kRef)
                kRef = kRank;
            else
                kRef = min(obj.kRef, nChSel);
            end

            % ── V_ref: calibration reference subspace ────────────────
            % Computed once from raw.calibration (channel-selected).
            % Used by energyCompaction and calibrationDrift groups.
            % V_ref_dom: [nChSel x kRef] dominant eigenvectors (sorted
            % descending by eigenvalue, matching Ux/Uy convention below).
            hasCalibration = isfield(obj.raw, 'calibration') && ...
                             ~isempty(obj.raw.calibration);
            if hasCalibration
                X_cal = double(obj.raw.calibration(chList, :));
                Cref  = TemporalAnalysis.robustCov( ...
                            X_cal, obj.useRobustCov, obj.blockSize);
                [V_ref_all, D_ref] = eig(Cref, 'vector');
                [~, ix_ref]        = sort(D_ref, 'descend');
                V_ref_dom          = V_ref_all(:, ix_ref(1:kRef));
                if obj.useRobustCov && ~isempty(which('block_geometric_median'))
                    covLabel = ' (robust cov)';
                else
                    covLabel = ' (standard cov)';
                end
                fprintf('  V_ref: calibration subspace [%d ch x %d dims]%s\n', ...
                    nChSel, kRef, covLabel);
            else
                V_ref_dom = [];
                fprintf('  V_ref: calibration segment absent — skipping calibrationDrift\n');
            end

            for c = 1:obj.nCombos
                fprintf('  combo %d/%d\n', c, obj.nCombos);

                % Stream combo on demand if a loader is wired up
                if ~isempty(obj.loaderFcn)
                    obj.loaderFcn(c);
                end

                sr = obj.subject_results(c);

                % Skip combos whose signals were not loaded
                if isempty(sr.cleanClosed) || isempty(sr.cleanOpen)
                    fprintf('  [skip] combo %d — signals not loaded\n', c);
                    continue;
                end

                cleanMap = struct('closed', double(sr.cleanClosed), ...
                                  'open',   double(sr.cleanOpen));

                segResults = struct();

                for si = 1:numel(obj.segments)
                    seg = obj.segments{si};

                    Xraw   = rawMap.(seg)(chList, :);
                    Xclean = cleanMap.(seg)(chList, :);

                    T = min(size(Xraw, 2), size(Xclean, 2));
                    Xraw   = Xraw(:,   1:T);
                    Xclean = Xclean(:, 1:T);

                    % Window start indices
                    starts = 1 : hopSamp : (T - winSamp + 1);
                    nWin   = numel(starts);

                    % Time axis: centre of each window in seconds
                    tAxis = (starts + winSamp/2 - 1) / obj.fs;   % [1 x nWin]

                    % Preallocate signal fields
                    sig = struct( ...
                        'rrmse',      zeros(1, nWin), ...
                        'corr',       zeros(1, nWin), ...
                        'snrDB',      zeros(1, nWin), ...
                        'rms_raw',    zeros(1, nWin), ...
                        'rms_clean',  zeros(1, nWin), ...
                        'std_raw',    zeros(1, nWin), ...
                        'std_clean',  zeros(1, nWin), ...
                        'mean_raw',   zeros(1, nWin), ...
                        'mean_clean', zeros(1, nWin));

                    % Preallocate subspace fields
                    sub = struct( ...
                        'angle_mean',     zeros(1, nWin), ...
                        'angle_max',      zeros(1, nWin), ...
                        'angle_all',      {cell(1, nWin)},  ...
                        'energy_ratio',   zeros(1, nWin), ...
                        'energy_removed', zeros(1, nWin), ...
                        'cov_fro_norm',   zeros(1, nWin), ...
                        'eigcorr',        zeros(1, nWin));

                    % Preallocate energy compaction fields
                    % (compute_energy_compaction_time.m)
                    ecp = struct( ...
                        'kCompact_raw',     zeros(1, nWin), ...
                        'kCompact_clean',   zeros(1, nWin), ...
                        'dimReduction',     zeros(1, nWin), ...
                        'energyCurve_raw',  {cell(1, nWin)}, ...
                        'energyCurve_clean',{cell(1, nWin)});

                    % Preallocate calibration drift fields
                    % (compute_energy_capture_time.m + compute_subspace_angle_time.m)
                    cal = struct( ...
                        'energyCapture_raw',   NaN(1, nWin), ...
                        'energyCapture_clean', NaN(1, nWin), ...
                        'energyCaptureGain',   NaN(1, nWin), ...
                        'angleFromRef_raw',    NaN(1, nWin), ...
                        'angleFromRef_clean',  NaN(1, nWin), ...
                        'angleFromRefGain',    NaN(1, nWin));



                    % ── Main window loop ─────────────────────────────
                    for w = 1:nWin
                        ix = starts(w) : starts(w) + winSamp - 1;

                        Xr = Xraw(:,   ix);
                        Xc = Xclean(:, ix);

                        % ── Signal stats ────────────────────────────
                        diff_w = Xc - Xr;

                        % RRMSE: rms(diff)/rms(raw) per channel, then mean
                        rms_r_ch = sqrt(mean(Xr.^2, 2));
                        rms_d_ch = sqrt(mean(diff_w.^2, 2));
                        rrmse_ch = rms_d_ch ./ max(rms_r_ch, eps);
                        sig.rrmse(w)     = mean(rrmse_ch, 'omitnan');

                        sig.rms_raw(w)   = mean(rms_r_ch, 'omitnan');
                        sig.rms_clean(w) = mean(sqrt(mean(Xc.^2, 2)), 'omitnan');
                        sig.std_raw(w)   = mean(std(Xr, 0, 2), 'omitnan');
                        sig.std_clean(w) = mean(std(Xc, 0, 2), 'omitnan');
                        sig.mean_raw(w)  = mean(Xr(:), 'omitnan');
                        sig.mean_clean(w)= mean(Xc(:), 'omitnan');

                        % SNR per channel, then mean
                        pow_r  = mean(Xr.^2, 2);
                        pow_d  = mean(diff_w.^2, 2);
                        snr_ch = 10 * log10(pow_r ./ max(pow_d, eps));
                        sig.snrDB(w) = mean(snr_ch, 'omitnan');

                        % Per-channel correlation, then mean
                        corr_ch = NaN(nChSel, 1);
                        for ch = 1:nChSel
                            xr = Xr(ch, :);
                            xc = Xc(ch, :);
                            if std(xr) > 0 && std(xc) > 0
                                rv = corrcoef(xr, xc);
                                corr_ch(ch) = rv(1,2);
                            end
                        end
                        sig.corr(w) = mean(corr_ch, 'omitnan');

                        % ── Subspace stats ───────────────────────────
                        % Remove per-channel mean before covariance
                        Xr0 = Xr - mean(Xr, 2);
                        Xc0 = Xc - mean(Xc, 2);

                        Cx = (Xr0 * Xr0') / (winSamp - 1);
                        Cy = (Xc0 * Xc0') / (winSamp - 1);

                        % Eigendecomposition — sorted descending
                        [Ux, Dx] = eig(Cx, 'vector');
                        [Uy, Dy] = eig(Cy, 'vector');
                        [Dx, ix_] = sort(Dx, 'descend');
                        [Dy, iy_] = sort(Dy, 'descend');
                        Ux = Ux(:, ix_);
                        Uy = Uy(:, iy_);

                        % Principal angles between top-k subspaces
                        ang = TemporalAnalysis.principalAngles(Ux(:,1:kRank), Uy(:,1:kRank));
                        sub.angle_mean(w) = mean(ang, 'omitnan');
                        sub.angle_max(w)  = max(ang);
                        sub.angle_all{w}  = ang;

                        % Energy ratio
                        sumDx = sum(Dx);
                        sumDy = sum(Dy);
                        if sumDx > 0
                            sub.energy_ratio(w)   = sumDy / sumDx;
                        else
                            sub.energy_ratio(w)   = NaN;
                        end
                        sub.energy_removed(w) = 1 - sub.energy_ratio(w);

                        % Covariance Frobenius distortion
                        normCx = norm(Cx, 'fro');
                        if normCx > 0
                            sub.cov_fro_norm(w) = norm(Cy - Cx, 'fro') / normCx;
                        else
                            sub.cov_fro_norm(w) = NaN;
                        end

                        % Eigenvalue spectrum correlation
                        m = min(numel(Dx), numel(Dy));
                        if m >= 2
                            ev = corrcoef(Dx(1:m), Dy(1:m));
                            sub.eigcorr(w) = ev(1,2);
                        else
                            sub.eigcorr(w) = NaN;
                        end

                        % ── BLOCK A: Energy compaction ───────────────
                        % From compute_energy_compaction_time.m
                        % Reuses Dx, Dy (already sorted descending).
                        % Dx, Dy are the eigenvalue vectors of Cx, Cy.
                        % Compute cumulative energy and effective dimensionality
                        % for raw and cleaned signals independently.

                        sumDx_safe = max(sumDx, eps);
                        sumDy_safe = max(sum(Dy), eps);

                        cumE_raw   = cumsum(Dx) / sumDx_safe;
                        cumE_clean = cumsum(Dy) / sumDy_safe;

                        k_raw = find(cumE_raw   >= obj.compactionThresh, 1);
                        k_cln = find(cumE_clean >= obj.compactionThresh, 1);

                        if isempty(k_raw),   k_raw = nChSel; end
                        if isempty(k_cln),   k_cln = nChSel; end

                        ecp.kCompact_raw(w)      = k_raw;
                        ecp.kCompact_clean(w)    = k_cln;
                        ecp.dimReduction(w)      = k_raw - k_cln;
                        ecp.energyCurve_raw{w}   = cumE_raw;
                        ecp.energyCurve_clean{w} = cumE_clean;

                        % ── BLOCK B: Calibration drift ───────────────
                        % From compute_energy_capture_time.m and
                        % compute_subspace_angle_time.m
                        % Requires V_ref_dom (precomputed from calibration).
                        % Reuses Cx, Cy, Ux, Uy from the subspace block above.

                        if ~isempty(V_ref_dom)

                            % Energy capture: fraction of window energy
                            % that lies in the calibration reference subspace
                            % P = V_ref' * C * V_ref  (trace = projected energy)
                            trCx = max(trace(Cx), eps);
                            trCy = max(trace(Cy), eps);

                            P_raw   = V_ref_dom' * Cx * V_ref_dom;
                            P_clean = V_ref_dom' * Cy * V_ref_dom;

                            cal.energyCapture_raw(w)   = trace(P_raw)   / trCx;
                            cal.energyCapture_clean(w) = trace(P_clean) / trCy;
                            cal.energyCaptureGain(w)   = ...
                                cal.energyCapture_clean(w) - cal.energyCapture_raw(w);

                            % Subspace angle vs reference
                            % principalAngles reuses the QR-normalised path
                            % Ux(:,1:kRef), Uy(:,1:kRef) are already sorted descending
                            ang_raw   = TemporalAnalysis.principalAngles( ...
                                            V_ref_dom, Ux(:,1:kRef));
                            ang_clean = TemporalAnalysis.principalAngles( ...
                                            V_ref_dom, Uy(:,1:kRef));

                            cal.angleFromRef_raw(w)   = mean(ang_raw,   'omitnan');
                            cal.angleFromRef_clean(w) = mean(ang_clean, 'omitnan');
                            cal.angleFromRefGain(w)   = ...
                                cal.angleFromRef_raw(w) - cal.angleFromRef_clean(w);
                        end

                    end % window loop

                    % ── Time-averaged summary scalars ────────────────
                    summary = struct();
                    summary.signal    = TemporalAnalysis.aggregateFields(sig);
                    summary.subspace  = TemporalAnalysis.aggregateFields( ...
                        rmfield(sub, 'angle_all'));   % skip cell field
                    summary.energyCompaction = TemporalAnalysis.aggregateFields( ...
                        rmfield(ecp, {'energyCurve_raw','energyCurve_clean'}));
                    summary.calibrationDrift = TemporalAnalysis.aggregateFields(cal);

                    segResults.(seg) = struct( ...
                        'time',             tAxis,  ...
                        'signal',           sig,    ...
                        'subspace',         sub,    ...
                        'energyCompaction', ecp,    ...
                        'calibrationDrift', cal,    ...
                        'summary',          summary);


                end % segment loop

                obj.results(c).segment    = segResults;
                obj.results(c).parameters = sr.parameters;
                obj.results(c).windowSec  = obj.windowSec;
                obj.results(c).hopSec     = obj.hopSec;
                obj.results(c).kRank      = kRank;
                obj.results(c).kRef       = kRef;
                obj.results(c).channels   = chList;

                % Free combo signals if a loader is managing memory
                if ~isempty(obj.unloaderFcn)
                    obj.unloaderFcn(c);
                end

            end % combo loop

            fprintf('  Done.\n');
        end

        % ── getSummary ───────────────────────────────────────────────
        function vals = getSummary(obj, fieldName, segment)
            % GETSUMMARY - Extract one time-averaged scalar across all combos.
            %
            %   fieldName : 'group.field_aggregation'
            %     group        = 'signal' | 'subspace'
            %     field        = any time-series field name
            %     aggregation  = 'mean' | 'std' | 'max' | 'p95'
            %
            %   Examples:
            %     ta.getSummary('signal.rrmse_mean',          'closed')
            %     ta.getSummary('signal.corr_mean',           'closed')
            %     ta.getSummary('signal.snrDB_mean',          'closed')
            %     ta.getSummary('signal.rrmse_p95',           'closed')
            %     ta.getSummary('subspace.angle_mean_mean',   'closed')
            %     ta.getSummary('subspace.angle_max_max',     'closed')
            %     ta.getSummary('subspace.cov_fro_norm_mean', 'closed')
            %     ta.getSummary('subspace.energy_removed_max','closed')
            %     ta.getSummary('subspace.eigcorr_mean',      'closed')
            %
            %   Returns [1 x nCombos] double — NaN for missing entries.

            obj.checkComputed();
            if nargin < 3, segment = 'closed'; end

            vals  = NaN(1, obj.nCombos);
            parts = strsplit(fieldName, '.');   % {'signal','rrmse_mean'} or {'subspace','cov_fro_norm_mean'}

            if numel(parts) ~= 2
                error('TemporalAnalysis:badField', ...
                    'fieldName must be ''group.field_aggregation'', e.g. ''signal.rrmse_mean''.');
            end

            grp   = parts{1};   % 'signal' or 'subspace'
            fkey  = parts{2};   % e.g. 'rrmse_mean' or 'cov_fro_norm_mean'

            for c = 1:obj.nCombos
                try
                    s = obj.results(c).segment.(segment).summary.(grp);
                    vals(c) = s.(fkey);
                catch
                    % field absent — leave NaN
                end
            end
        end

        % ── getTimeSeries ────────────────────────────────────────────
        function [t, v] = getTimeSeries(obj, comboIdx, fieldName, segment)
            % GETTIMESERIES - Return time axis and one windowed series for a combo.
            %
            %   comboIdx  : integer in [1, nCombos]
            %   fieldName : 'group.field'  e.g. 'signal.rrmse', 'subspace.angle_mean'
            %   segment   : 'closed' (default) | 'open'
            %
            %   Returns:
            %     t [1 x W]  window centre times (sec)
            %     v [1 x W]  values

            obj.checkComputed();
            if nargin < 4, segment = 'closed'; end

            parts = strsplit(fieldName, '.');
            if numel(parts) ~= 2
                error('TemporalAnalysis:badField', ...
                    'fieldName must be ''group.field'', e.g. ''signal.rrmse''.');
            end

            seg_r = obj.results(comboIdx).segment.(segment);
            t = seg_r.time;
            v = seg_r.(parts{1}).(parts{2});
        end

        % ── paramTable ───────────────────────────────────────────────
        function T = paramTable(obj, segment)
            % PARAMTABLE - One row per combo: parameters + key temporal scalars.
            obj.checkComputed();
            if nargin < 2, segment = 'closed'; end

            rows = cell(obj.nCombos, 1);
            for c = 1:obj.nCombos
                row = table();
                row.comboIdx = c;

                params = obj.results(c).parameters;
                pf = fieldnames(params);
                for i = 1:numel(pf)
                    val = params.(pf{i});
                    if isscalar(val) && isnumeric(val)
                        row.(pf{i}) = val;
                    end
                end

                if isfield(obj.results(c).segment, segment)
                    s = obj.results(c).segment.(segment).summary;

                    % Signal group
                    row.rrmse_mean  = s.signal.rrmse_mean;
                    row.rrmse_p95   = s.signal.rrmse_p95;
                    row.corr_mean   = s.signal.corr_mean;
                    row.snrDB_mean  = s.signal.snrDB_mean;

                    % Subspace group
                    row.angle_mean_mean    = s.subspace.angle_mean_mean;
                    row.angle_max_max      = s.subspace.angle_max_max;
                    row.cov_fro_mean       = s.subspace.cov_fro_norm_mean;
                    row.energy_removed_mean= s.subspace.energy_removed_mean;
                    row.eigcorr_mean       = s.subspace.eigcorr_mean;

                    % Energy compaction group
                    if isfield(s, 'energyCompaction')
                        ec = s.energyCompaction;
                        row.kCompact_raw_mean   = ec.kCompact_raw_mean;
                        row.kCompact_clean_mean = ec.kCompact_clean_mean;
                        row.dimReduction_mean   = ec.dimReduction_mean;
                        row.dimReduction_p95    = ec.dimReduction_p95;
                    end

                    % Calibration drift group
                    if isfield(s, 'calibrationDrift')
                        cd = s.calibrationDrift;
                        row.energyCaptureGain_mean  = cd.energyCaptureGain_mean;
                        row.energyCapture_raw_mean  = cd.energyCapture_raw_mean;
                        row.energyCapture_clean_mean= cd.energyCapture_clean_mean;
                        row.angleFromRefGain_mean   = cd.angleFromRefGain_mean;
                        row.angleFromRef_raw_mean   = cd.angleFromRef_raw_mean;
                        row.angleFromRef_clean_mean = cd.angleFromRef_clean_mean;
                    end
                end

                rows{c} = row;
            end

            T = vertcat(rows{:});
        end

        % ── plotTimeSeries ───────────────────────────────────────────
        function plotTimeSeries(obj, comboIdx, fieldName, segment)
            % PLOTTIMESERIES - Plot one windowed series for a single combo.
            %
            %   ta.plotTimeSeries(1, 'signal.rrmse',       'closed')
            %   ta.plotTimeSeries(1, 'subspace.angle_mean','open')

            obj.checkComputed();
            if nargin < 3, fieldName = 'signal.rrmse'; end
            if nargin < 4, segment   = 'closed';        end

            [t, v] = obj.getTimeSeries(comboIdx, fieldName, segment);

            figure('Name', sprintf('%s — %s (S%d combo %d %s)', ...
                obj.algorithmName, fieldName, obj.subjectID, comboIdx, segment));
            plot(t, v, 'LineWidth', 1.4);
            xlabel('Time (s)'); ylabel(strrep(fieldName, '.', ' — '));
            title(sprintf('%s | S%d combo %d | %s', ...
                strrep(fieldName,'_',' '), obj.subjectID, comboIdx, segment), ...
                'Interpreter','none');
            grid on;
        end

        % ── plotComparison ───────────────────────────────────────────
        function plotComparison(obj, fieldName, segment)
            % PLOTCOMPARISON - Overlay all combos for one field (parameter sweep view).
            %
            %   ta.plotComparison('signal.rrmse',        'closed')
            %   ta.plotComparison('subspace.angle_mean', 'closed')

            obj.checkComputed();
            if nargin < 2, fieldName = 'signal.rrmse'; end
            if nargin < 3, segment   = 'closed';        end

            figure('Name', sprintf('%s — %s %s all combos', ...
                obj.algorithmName, fieldName, segment));
            hold on;
            legLabels = cell(obj.nCombos, 1);

            for c = 1:obj.nCombos
                [t, v] = obj.getTimeSeries(c, fieldName, segment);
                plot(t, v, 'LineWidth', 1.2);
                legLabels{c} = sprintf('combo %d', c);
            end

            xlabel('Time (s)');
            ylabel(strrep(fieldName, '.', ' — '));
            title(sprintf('%s | S%d | %s | %s — all combos', ...
                obj.algorithmName, obj.subjectID, segment, ...
                strrep(fieldName,'_',' ')), 'Interpreter','none');
            legend(legLabels, 'Location','northeastoutside');
            grid on;
        end

        % ── plotSignalPanel ──────────────────────────────────────────
        function plotSignalPanel(obj, comboIdx, segment)
            % PLOTSIGNALPANEL - 3-panel summary: rrmse, corr, snrDB vs time.

            obj.checkComputed();
            if nargin < 2, comboIdx = 1;       end
            if nargin < 3, segment  = 'closed'; end

            seg_r = obj.results(comboIdx).segment.(segment);
            t     = seg_r.time;
            sig   = seg_r.signal;

            figure('Name', sprintf('%s — signal panel S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));

            subplot(3,1,1);
            plot(t, sig.rrmse, 'LineWidth', 1.3); ylabel('RRMSE'); grid on;
            title(sprintf('%s — S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter','none');
            yline(0, 'k--');

            subplot(3,1,2);
            plot(t, sig.corr,  'LineWidth', 1.3); ylabel('Correlation'); grid on;
            yline(1, 'k--'); ylim([max(-1, min(sig.corr)-0.05), 1.05]);

            subplot(3,1,3);
            plot(t, sig.snrDB, 'LineWidth', 1.3); ylabel('SNR (dB)'); grid on;
            xlabel('Time (s)');
        end

        % ── plotSubspacePanel ────────────────────────────────────────
        function plotSubspacePanel(obj, comboIdx, segment)
            % PLOTSUBSPACEPANEL - 4-panel subspace summary vs time.

            obj.checkComputed();
            if nargin < 2, comboIdx = 1;       end
            if nargin < 3, segment  = 'closed'; end

            seg_r = obj.results(comboIdx).segment.(segment);
            t     = seg_r.time;
            sub   = seg_r.subspace;

            figure('Name', sprintf('%s — subspace panel S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));

            subplot(4,1,1);
            plot(t, sub.angle_mean, 'LineWidth', 1.3);
            hold on;
            plot(t, sub.angle_max,  '--', 'LineWidth', 1);
            legend({'mean angle','max angle'},'Location','northeast');
            ylabel('Principal angle (°)'); grid on;
            title(sprintf('%s — subspace | S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter','none');

            subplot(4,1,2);
            plot(t, sub.energy_removed, 'LineWidth', 1.3);
            ylabel('Energy removed'); grid on;
            yline(0,'k--');

            subplot(4,1,3);
            plot(t, sub.cov_fro_norm, 'LineWidth', 1.3);
            ylabel('Cov Frobenius norm'); grid on;
            yline(0,'k--');

            subplot(4,1,4);
            plot(t, sub.eigcorr, 'LineWidth', 1.3);
            ylabel('Eigenvalue corr'); grid on;
            xlabel('Time (s)');
            yline(1,'k--'); ylim([max(-1, min(sub.eigcorr)-0.05), 1.05]);
        end

        % ── plotCompactionPanel ───────────────────────────────────────
        function plotCompactionPanel(obj, comboIdx, segment)
            % PLOTCOMPACTIONPANEL - 3-panel energy compaction summary vs time.
            %
            %   Shows how many dimensions are needed to capture compactionThresh%
            %   of signal energy for raw vs cleaned, and the reduction.
            %
            %   ta.plotCompactionPanel(1, 'closed')

            obj.checkComputed();
            if nargin < 2, comboIdx = 1;       end
            if nargin < 3, segment  = 'closed'; end

            seg_r = obj.results(comboIdx).segment.(segment);
            t     = seg_r.time;
            ecp   = seg_r.energyCompaction;

            figure('Name', sprintf('%s — energy compaction S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));

            subplot(3,1,1);
            plot(t, ecp.kCompact_raw,   'LineWidth', 1.3); hold on;
            plot(t, ecp.kCompact_clean, '--', 'LineWidth', 1.3);
            legend({'raw','clean'}, 'Location','northeast');
            ylabel(sprintf('Dims for %.0f%% energy', obj.compactionThresh*100));
            title(sprintf('%s — energy compaction | S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter','none');
            grid on;

            subplot(3,1,2);
            bar(t, ecp.dimReduction, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor','none');
            ylabel('Dim reduction (raw - clean)');
            yline(0,'k--');
            grid on;

            subplot(3,1,3);
            % Show mean compaction curve for first and last window (illustrative)
            nW = numel(ecp.energyCurve_raw);
            if nW > 0
                nC = numel(ecp.energyCurve_raw{1});
                xax = 1:nC;
                plot(xax, ecp.energyCurve_raw{1},   'b-',  'LineWidth', 1.2); hold on;
                plot(xax, ecp.energyCurve_clean{1}, 'b--', 'LineWidth', 1.2);
                if nW > 1
                    plot(xax, ecp.energyCurve_raw{end},   'r-',  'LineWidth', 1.2);
                    plot(xax, ecp.energyCurve_clean{end}, 'r--', 'LineWidth', 1.2);
                    legend({'raw w1','clean w1','raw wN','clean wN'},'Location','southeast');
                else
                    legend({'raw','clean'},'Location','southeast');
                end
                yline(obj.compactionThresh, 'k:', sprintf('%.0f%%', obj.compactionThresh*100));
            end
            xlabel('Component index'); ylabel('Cumulative energy');
            grid on;
        end

        % ── plotCalibrationDrift ─────────────────────────────────────
        function plotCalibrationDrift(obj, comboIdx, segment)
            % PLOTCALIBRATIONDRIFT - 4-panel calibration drift summary vs time.
            %
            %   Shows how much the signal departs from the calibration reference
            %   subspace window by window, and whether ASR corrects the drift.
            %
            %     Panel 1: energyCapture_raw vs clean (fraction in calib subspace)
            %     Panel 2: energyCaptureGain  (>0 = ASR pushed toward calibration)
            %     Panel 3: angleFromRef raw vs clean (deg from calib subspace)
            %     Panel 4: angleFromRefGain   (>0 = ASR moved closer to calibration)
            %
            %   ta.plotCalibrationDrift(1, 'closed')

            obj.checkComputed();
            if nargin < 2, comboIdx = 1;       end
            if nargin < 3, segment  = 'closed'; end

            seg_r = obj.results(comboIdx).segment.(segment);
            t     = seg_r.time;
            cal   = seg_r.calibrationDrift;

            if all(isnan(cal.energyCapture_raw))
                warning('TemporalAnalysis:noCalibration', ...
                    'calibrationDrift is all NaN (no calibration segment was found).');
                return;
            end

            figure('Name', sprintf('%s — calibration drift S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));

            subplot(4,1,1);
            plot(t, cal.energyCapture_raw,   'LineWidth', 1.3); hold on;
            plot(t, cal.energyCapture_clean, '--', 'LineWidth', 1.3);
            legend({'raw','clean'}, 'Location','best');
            ylabel('Energy in calib subspace');
            title(sprintf('%s — calibration drift | S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter','none');
            ylim([0 1.05]); grid on;

            subplot(4,1,2);
            bar(t, cal.energyCaptureGain, ...
                'FaceColor', [0.2 0.7 0.3], 'EdgeColor','none');
            ylabel('Energy capture gain (clean-raw)');
            yline(0,'k--'); grid on;

            subplot(4,1,3);
            plot(t, cal.angleFromRef_raw,   'LineWidth', 1.3); hold on;
            plot(t, cal.angleFromRef_clean, '--', 'LineWidth', 1.3);
            legend({'raw','clean'}, 'Location','best');
            ylabel('Angle from calib ref (°)'); grid on;

            subplot(4,1,4);
            bar(t, cal.angleFromRefGain, ...
                'FaceColor', [0.8 0.4 0.2], 'EdgeColor','none');
            ylabel('Angle gain raw-clean (°)');
            yline(0,'k--'); grid on;
            xlabel('Time (s)');
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function checkComputed(obj)
            if isempty(obj.results) || ~isfield(obj.results(1), 'segment')
                error('TemporalAnalysis:notComputed', ...
                    'Call compute() before accessing results.');
            end
        end

    end % private methods

    % ================================================================
    methods (Static, Access = public)

        % ── principalAngles ──────────────────────────────────────────
        function ang = principalAngles(U, V)
            % PRINCIPALANGLES - Principal angles between column spaces of U and V.
            %   Returns angles in degrees, sorted ascending.
            %   Each column of U and V should be a unit vector (orthonormal basis).
            %
            %   Implementation: SVD of the cross-Gram matrix U'*V.
            %   The singular values are cos(theta_i) for principal angle theta_i.
            %
            %   Note: eig() eigenvectors are NOT guaranteed to be orthonormal when
            %   eigenvalues are repeated.  We orthonormalise via QR before computing
            %   angles so the result is numerically stable.

            if isempty(U) || isempty(V)
                ang = NaN;
                return;
            end

            % Orthonormalize via economy QR — makes SVD numerically meaningful
            [Uq, ~] = qr(U, 'econ');
            [Vq, ~] = qr(V, 'econ');

            s = svd(Uq' * Vq);

            % Clamp to [-1, 1] against floating point overshoot
            s = min(max(s, -1), 1);

            ang = sort(acosd(s), 'ascend');
        end

        % ── aggregateFields ──────────────────────────────────────────
        function out = aggregateFields(s)
            % AGGREGATEFIELDS - Compute mean/std/max/p95 for each numeric vector
            %   field in struct s.  Returns a flat struct with suffix-named fields.
            %
            %   Input:  s.(field) = [1 x W] numeric vector
            %   Output: out.(field_mean), out.(field_std), out.(field_max),
            %           out.(field_p95)

            out   = struct();
            flds  = fieldnames(s);

            for i = 1:numel(flds)
                f = flds{i};
                v = s.(f);

                if ~isnumeric(v) || ~isvector(v)
                    continue;   % skip cell arrays (angle_all) and non-vectors
                end

                v = double(v(:));
                vv = v(~isnan(v));   % strip NaNs once — std/max/prctile don't all support 'omitnan'
                out.([f '_mean']) = mean(v,  'omitnan');
                out.([f '_std'])  = std(vv);
                out.([f '_max'])  = max(vv);
                out.([f '_p95'])  = prctile(vv, 95);
            end
        end

        % ── robustCov ────────────────────────────────────────────────
        function C = robustCov(X, useRobust, blockSize)
            % ROBUSTCOV - Robust or standard covariance of [nCh x T] signal.
            %
            %   Implements the block geometric median approach from
            %   compute_state_covariance_rms.m when block_geometric_median
            %   (EEGLAB) is available.  Falls back to standard cov() silently.
            %
            %   X         : [C x T]  EEG data
            %   useRobust : logical  (true = try block_geometric_median)
            %   blockSize : integer  subsampling factor (default 10)
            %
            %   Returns C: [C x C] covariance matrix

            if nargin < 2, useRobust = true; end
            if nargin < 3, blockSize = 10;   end

            [Ch, T] = size(X);

            if useRobust && ~isempty(which('block_geometric_median'))
                try
                    % Accumulate block outer products
                    % (mirrors compute_state_covariance_rms.m exactly)
                    Xt = X';   % [T x C]
                    nBlocks = length(1:blockSize:T);
                    U = zeros(nBlocks, Ch*Ch);

                    for k = 1:blockSize
                        range = min(T, k:blockSize:(T+k-1));
                        U = U + reshape( ...
                            bsxfun(@times, ...
                                reshape(Xt(range,:), [], 1, Ch), ...
                                reshape(Xt(range,:), [], Ch, 1)), ...
                            size(U));
                    end

                    C = real(reshape(block_geometric_median(U / blockSize), Ch, Ch));
                    return;
                catch
                    % fall through to standard cov
                end
            end

            % Standard covariance (MATLAB cov normalises by T-1)
            C = cov(X');
        end

    end % static methods

end
