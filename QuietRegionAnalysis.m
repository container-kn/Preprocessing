classdef QuietRegionAnalysis < handle
    % QuietRegionAnalysis - Signal preservation metrics on blink-free regions.
    %
    %   Whole-signal metrics (SignalStatistics, SpectralAnalysis) conflate
    %   artefact-contaminated samples with clean background EEG.  This class
    %   isolates the NON-blink samples and evaluates how faithfully ASR
    %   preserved the underlying signal where nothing needed correcting.
    %
    %   ── Quiet mask construction ──────────────────────────────────────────
    %   For each channel, every raw blink index is expanded by ±quietWinSec
    %   (default 0.25 s) to form an exclusion zone.  Samples outside ALL
    %   exclusion zones are "quiet".  The mask is built from raw blink indices
    %   so it is the same across all parameter combos (fair comparison).
    %
    %   Requires a computed BlinkAnalysis instance — uses its per-channel
    %   idx_raw_closed / idx_raw_open cell arrays.
    %
    %   ── Per-channel metrics on quiet samples ─────────────────────────────
    %     corr         Pearson correlation between orig_quiet and clean_quiet
    %     rre          Relative RMS error: rms(clean-orig)/rms(orig)  [= RRMSE]
    %     varRatio     var(clean) / var(orig)  — variance preservation
    %     residualDelay  finddelay(orig, clean, ±maxDelaySec)  in samples
    %                    positive = clean is delayed; 0 = perfect alignment
    %
    %   ── Band power change on quiet regions ───────────────────────────────
    %   SpectralAnalysis computes band power on the full signal (contaminated).
    %   Here we compute it on quiet samples only — the meaningful version.
    %     bandDelta.<band>  (pClean - pOrig) / pOrig  per channel
    %
    %   Bands: delta[1-4], theta[4-8], alpha[8-13], beta[13-30], gamma[30-100] Hz
    %
    %   ── Summary scalars (results(k).summary.(segment)) ──────────────────
    %     corr_median / corr_mean / corr_per_channel
    %     rre_median  / rre_mean  / rre_per_channel
    %     varRatio_median / varRatio_mean / varRatio_per_channel
    %     residualDelay_median_samples / residualDelay_per_channel
    %     bandDelta_<band>_median / bandDelta_<band>_per_channel
    %     quietFraction     fraction of samples that were quiet
    %     nQuietSamples     absolute count of quiet samples used
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       bl = BlinkAnalysis(...); bl.compute();
    %       qr = QuietRegionAnalysis(raw, subject_results, fs, bl, subjectID, algo);
    %       qr.compute();
    %       qr.getSummary('corr_median',    'closed')
    %       qr.getSummary('rre_median',     'closed')
    %       qr.getSummary('varRatio_median','closed')
    %       T = qr.paramTable();
    %
    %   ── Relationship to computeASRMetrics.m ─────────────────────────────
    %   Faithful OOP rewrite of computeASRMetrics + local_bandpower:
    %     corr      → M.quiet.corr_median
    %     rre       → M.quiet.rre_median
    %     varRatio  → M.quiet.varRatio_median
    %     delay     → M.quiet.residualDelay_median_samples
    %     bpDelta   → M.bandChanges.relChange_median
    %   Added: per-channel mean/median/std, band power on quiet (not whole signal),
    %          quietFraction, nQuietSamples, open-eyes segment.

    % ================================================================
    properties (Access = public)

        raw               % struct: .calibration/.closed/.open [C x T]
        subject_results   % 1xN struct array (.cleanClosed/.cleanOpen)
        fs
        subjectID
        algorithmName
        nCombos

        blinkAna          % BlinkAnalysis instance (must be computed)

        % Config
        quietWinSec  = 0.25   % ± exclusion window around each blink (sec)
        maxDelaySec  = 0.25   % finddelay search limit (sec)
        minQuietSamples = 50  % skip channel if fewer quiet samples than this
        segments = {'closed', 'open'}

        bands = struct( ...
            'delta', [1  4],  ...
            'theta', [4  8],  ...
            'alpha', [8  13], ...
            'beta',  [13 30], ...
            'gamma', [30 100])

        % Results — 1xN struct array
        %   results(k).summary.(segment)          — scalar summaries
        %   results(k).perChannel.(segment)       — per-channel arrays
        %   results(k).quietMask.(segment)        — [1 x T] logical mask (shared raw)
        %   results(k).parameters
        results

        % Optional streaming loader — set by ExperimentAnalysis for split-file layout.
        loaderFcn   = []   % @(k) ana.loadCombo(k)
        unloaderFcn = []   % @(k) ana.unloadCombo(k)

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = QuietRegionAnalysis(raw, subject_results, fs, blinkAna, subjectID, algorithmName)
            % QUIETREGIONANALYSIS
            %   raw            : struct with .calibration/.closed/.open [C x T]
            %   subject_results: 1xN struct array (.cleanClosed/.cleanOpen)
            %   fs             : sampling rate Hz
            %   blinkAna       : BlinkAnalysis instance (must have been compute()d)
            %   subjectID      : scalar (default 0)
            %   algorithmName  : string (default 'unknown')

            if nargin < 5, subjectID     = 0;         end
            if nargin < 6, algorithmName = 'unknown'; end

            if ~isa(blinkAna, 'BlinkAnalysis')
                error('QuietRegionAnalysis:badInput', ...
                    'blinkAna must be a BlinkAnalysis instance.');
            end

            obj.raw             = raw;
            obj.subject_results = subject_results;
            obj.fs              = fs;
            obj.blinkAna        = blinkAna;
            obj.subjectID       = subjectID;
            obj.algorithmName   = algorithmName;
            obj.nCombos         = numel(subject_results);
            obj.results         = struct();
        end

        % ── compute ──────────────────────────────────────────────────
        function compute(obj, varargin)
            % COMPUTE - Run quiet-region analysis across all combos.
            %
            %   Optional name-value:
            %     'quietWinSec'      0.25   exclusion half-window around blinks (sec)
            %     'maxDelaySec'      0.25   finddelay search range (sec)
            %     'minQuietSamples'  50     skip channel if fewer quiet samples

            p = inputParser;
            addParameter(p, 'quietWinSec',     obj.quietWinSec,     @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'maxDelaySec',     obj.maxDelaySec,     @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'minQuietSamples', obj.minQuietSamples, @(x) isnumeric(x) && isscalar(x));
            parse(p, varargin{:});
            obj.quietWinSec     = p.Results.quietWinSec;
            obj.maxDelaySec     = p.Results.maxDelaySec;
            obj.minQuietSamples = p.Results.minQuietSamples;

            obj.checkBlinkComputed();

            fprintf('QuietRegionAnalysis: S%d %s (%d combo(s))...\n', ...
                obj.subjectID, obj.algorithmName, obj.nCombos);

            Q       = max(1, round(obj.quietWinSec * obj.fs));
            maxLag  = round(obj.maxDelaySec * obj.fs);

            % Raw signals per segment (same across all combos)
            rawMap = struct('closed', double(obj.raw.closed), ...
                            'open',   double(obj.raw.open));

            for c = 1:obj.nCombos
                % Stream combo on demand if a loader is wired up
                if ~isempty(obj.loaderFcn)
                    obj.loaderFcn(c);
                end

                sr = obj.subject_results(c);

                if isempty(sr.cleanClosed) || isempty(sr.cleanOpen)
                    fprintf('  [skip] combo %d — signals not loaded\n', c);
                    continue;
                end

                cleanMap = struct('closed', double(sr.cleanClosed), ...
                                  'open',   double(sr.cleanOpen));

                perCh   = struct();
                summ    = struct();
                qmasks  = struct();

                for si = 1:numel(obj.segments)
                    seg = obj.segments{si};

                    Xraw  = rawMap.(seg);
                    Xclean = cleanMap.(seg);

                    T  = min(size(Xraw,2), size(Xclean,2));
                    Xraw   = Xraw(:, 1:T);
                    Xclean = Xclean(:, 1:T);
                    nCh    = size(Xraw, 1);

                    % Blink indices for this segment from BlinkAnalysis
                    % (these are from combo c — same raw mask is fair because
                    %  the threshold is raw-fixed, but we use combo c's indices
                    %  so the mask reflects what the detector found for this run)
                    idxBlinkField = sprintf('idx_raw_%s', seg);
                    blinkIdxCells = obj.blinkAna.results(c).perChannel.(idxBlinkField);

                    % Per-channel preallocate
                    corr_ch     = NaN(nCh, 1);
                    rre_ch      = NaN(nCh, 1);
                    varRatio_ch = NaN(nCh, 1);
                    delay_ch    = NaN(nCh, 1);
                    bnames      = fieldnames(obj.bands);
                    bdelta_ch   = NaN(nCh, numel(bnames));
                    nQuiet_ch   = zeros(nCh, 1);

                    % Shared quiet mask: build from ALL channels' blinks combined
                    % so that noisy channels don't inflate quiet time for clean ones.
                    % This matches computeASRMetrics.m approach (per-channel mask)
                    % but we also expose a shared mask for plotting convenience.
                    sharedMask = true(1, T);
                    for ch = 1:nCh
                        bIdx = blinkIdxCells{ch}(:);
                        for m = 1:numel(bIdx)
                            j = bIdx(m);
                            if j >= 1 && j <= T
                                sharedMask(max(1,j-Q) : min(T,j+Q)) = false;
                            end
                        end
                    end

                    for ch = 1:nCh
                        % Per-channel quiet mask from this channel's blinks only
                        % (matches computeASRMetrics.m exactly)
                        bIdx   = blinkIdxCells{ch}(:);
                        qMask  = true(1, T);
                        for m = 1:numel(bIdx)
                            j = bIdx(m);
                            if j >= 1 && j <= T
                                qMask(max(1,j-Q) : min(T,j+Q)) = false;
                            end
                        end

                        nQ = sum(qMask);
                        nQuiet_ch(ch) = nQ;

                        if nQ < obj.minQuietSamples
                            continue;   % not enough clean background — skip
                        end

                        xo = double(Xraw(ch,   qMask));
                        yc = double(Xclean(ch, qMask));

                        % Correlation
                        if std(xo) > 0 && std(yc) > 0
                            corr_ch(ch) = corr(xo', yc');
                        end

                        % Relative RMS error  (= RRMSE on quiet samples)
                        rmsOrig = rms(xo);
                        if rmsOrig > 0
                            rre_ch(ch) = rms(yc - xo) / rmsOrig;
                        end

                        % Variance ratio
                        vOrig = var(xo);
                        if vOrig > 0
                            varRatio_ch(ch) = var(yc) / vOrig;
                        end

                        % Residual delay via cross-correlation
                        % Positive = yc (clean) is delayed relative to xo (orig)
                        delay_ch(ch) = finddelay(xo, yc, maxLag);

                        % Band power delta on QUIET samples only
                        for bi = 1:numel(bnames)
                            pO = QuietRegionAnalysis.bandpowerSafe(xo, obj.fs, obj.bands.(bnames{bi}));
                            pC = QuietRegionAnalysis.bandpowerSafe(yc, obj.fs, obj.bands.(bnames{bi}));
                            if pO > 0
                                bdelta_ch(ch, bi) = (pC - pO) / pO;
                            end
                        end
                    end

                    % --- Per-channel struct ---
                    bdeltaStruct = struct();
                    for bi = 1:numel(bnames)
                        bdeltaStruct.(bnames{bi}) = bdelta_ch(:, bi);
                    end

                    perCh.(seg) = struct( ...
                        'corr',          corr_ch,     ...
                        'rre',           rre_ch,      ...
                        'varRatio',      varRatio_ch, ...
                        'residualDelay', delay_ch,    ...
                        'bandDelta',     bdeltaStruct, ...
                        'nQuietSamples', nQuiet_ch);

                    % --- Summary scalars (nan-aware) ---
                    s = struct();
                    s.corr_median       = median(corr_ch,    'omitnan');
                    s.corr_mean         = mean(corr_ch,      'omitnan');
                    s.corr_std          = std(corr_ch,       'omitnan');
                    s.rre_median        = median(rre_ch,     'omitnan');
                    s.rre_mean          = mean(rre_ch,       'omitnan');
                    s.rre_std           = std(rre_ch,        'omitnan');
                    s.varRatio_median   = median(varRatio_ch,'omitnan');
                    s.varRatio_mean     = mean(varRatio_ch,  'omitnan');
                    s.varRatio_std      = std(varRatio_ch,   'omitnan');
                    s.residualDelay_median_samples = median(delay_ch, 'omitnan');
                    s.residualDelay_median_ms      = median(delay_ch, 'omitnan') / obj.fs * 1000;
                    s.nQuietSamples     = sum(nQuiet_ch);
                    s.quietFraction     = sum(sharedMask) / T;

                    for bi = 1:numel(bnames)
                        b = bnames{bi};
                        col = bdelta_ch(:, bi);
                        s.(sprintf('bandDelta_%s_median', b)) = median(col, 'omitnan');
                        s.(sprintf('bandDelta_%s_mean',   b)) = mean(col,   'omitnan');
                    end

                    summ.(seg)   = s;
                    qmasks.(seg) = sharedMask;
                end

                obj.results(c).perChannel = perCh;
                obj.results(c).summary    = summ;
                obj.results(c).quietMask  = qmasks;
                obj.results(c).parameters = sr.parameters;

                % Free combo signals if a loader is managing memory
                if ~isempty(obj.unloaderFcn)
                    obj.unloaderFcn(c);
                end
            end

            fprintf('  Done.\n');
        end

        % ── getSummary ───────────────────────────────────────────────
        function vals = getSummary(obj, fieldName, segment)
            % GETSUMMARY - Extract one scalar across all combos.
            %   fieldName : e.g. 'corr_median', 'rre_median', 'varRatio_median',
            %               'residualDelay_median_samples', 'quietFraction',
            %               'bandDelta_alpha_median'
            %   segment   : 'closed' (default) | 'open'
            obj.checkComputed();
            if nargin < 3, segment = 'closed'; end

            vals = NaN(1, obj.nCombos);
            for c = 1:obj.nCombos
                if isfield(obj.results(c).summary, segment) && ...
                   isfield(obj.results(c).summary.(segment), fieldName)
                    vals(c) = obj.results(c).summary.(segment).(fieldName);
                end
            end
        end

        % ── paramTable ───────────────────────────────────────────────
        function T = paramTable(obj, segment)
            % PARAMTABLE - One row per combo: parameters + quiet-region scalars.
            obj.checkComputed();
            if nargin < 2, segment = 'closed'; end

            rows   = cell(obj.nCombos, 1);
            bnames = fieldnames(obj.bands);

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

                if isfield(obj.results(c).summary, segment)
                    s = obj.results(c).summary.(segment);
                    row.corr_median          = s.corr_median;
                    row.rre_median           = s.rre_median;
                    row.varRatio_median      = s.varRatio_median;
                    row.delay_median_samples = s.residualDelay_median_samples;
                    row.delay_median_ms      = s.residualDelay_median_ms;
                    row.quietFraction        = s.quietFraction;
                    row.nQuietSamples        = s.nQuietSamples;
                    for bi = 1:numel(bnames)
                        b = bnames{bi};
                        row.(sprintf('%s_delta_median', b)) = s.(sprintf('bandDelta_%s_median', b));
                    end
                end

                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

        % ── plotQuietMetrics ─────────────────────────────────────────
        function plotQuietMetrics(obj, segment)
            % PLOTQUIETMETRICS - Box-like bar chart of corr/rre/varRatio across combos.
            obj.checkComputed();
            if nargin < 2, segment = 'closed'; end

            nC    = obj.nCombos;
            corrs = obj.getSummary('corr_median',    segment);
            rres  = obj.getSummary('rre_median',     segment);
            vrs   = obj.getSummary('varRatio_median',segment);

            figure('Name', sprintf('%s — quiet region metrics (%s)', ...
                obj.algorithmName, segment));

            subplot(1,3,1);
            bar(corrs); ylabel('Correlation (median)'); xlabel('Combo');
            title('Signal correlation'); ylim([0 1]); grid on;
            yline(1,'k--'); set(gca,'XTick',1:nC);

            subplot(1,3,2);
            bar(rres); ylabel('RRE (median)'); xlabel('Combo');
            title('Relative RMS error'); grid on;
            yline(0,'k--'); set(gca,'XTick',1:nC);

            subplot(1,3,3);
            bar(vrs); ylabel('Var ratio (median)'); xlabel('Combo');
            title('Variance preservation'); grid on;
            yline(1,'k--'); set(gca,'XTick',1:nC);

            sgtitle(sprintf('%s — quiet-region preservation (S%d, %s)', ...
                obj.algorithmName, obj.subjectID, segment), 'Interpreter','none');
        end

        % ── plotBandDeltas ────────────────────────────────────────────
        function plotBandDeltas(obj, comboIdx, segment)
            % PLOTBANDDELTAS - Per-channel band power delta (quiet regions only).
            obj.checkComputed();
            if nargin < 2, comboIdx = 1; end
            if nargin < 3, segment  = 'closed'; end

            bnames = fieldnames(obj.bands);
            nB     = numel(bnames);
            pc     = obj.results(comboIdx).perChannel.(segment);

            data = zeros(size(pc.corr, 1), nB);
            for bi = 1:nB
                data(:,bi) = pc.bandDelta.(bnames{bi});
            end

            figure('Name', sprintf('%s — quiet band Δ S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));
            bar(data, 'grouped');
            set(gca, 'XTick', 1:size(data,1));
            xlabel('Channel'); ylabel('(P_{clean} - P_{orig}) / P_{orig}');
            legend(bnames, 'Location','northeastoutside');
            yline(0, 'k--', 'LineWidth', 1.2);
            title(sprintf('%s — band power change on quiet regions (S%d, combo %d, %s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), 'Interpreter','none');
            grid on;
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function checkComputed(obj)
            if isempty(obj.results) || ~isfield(obj.results(1), 'summary')
                error('QuietRegionAnalysis:notComputed', ...
                    'Call compute() before accessing results.');
            end
        end

        function checkBlinkComputed(obj)
            if isempty(obj.blinkAna) || isempty(obj.blinkAna.results) || ...
               ~isfield(obj.blinkAna.results(1), 'perChannel')
                error('QuietRegionAnalysis:blinkNotComputed', ...
                    'BlinkAnalysis must be computed before QuietRegionAnalysis.');
            end
            if obj.blinkAna.nCombos ~= obj.nCombos
                error('QuietRegionAnalysis:combomismatch', ...
                    'BlinkAnalysis has %d combos but subject_results has %d.', ...
                    obj.blinkAna.nCombos, obj.nCombos);
            end
        end

    end % private methods

    % ================================================================
    methods (Static, Access = private)

        function p = bandpowerSafe(x, fs, band)
            % BANDPOWERSAFE - Band power via pwelch; falls back gracefully.
            %   Matches local_bandpower in computeASRMetrics.m:
            %   uses MATLAB's bandpower() if available, else pwelch+trapz.

            if numel(x) < 4
                p = 0;
                return;
            end

            if exist('bandpower', 'file') == 2
                try
                    p = bandpower(double(x), fs, band);
                    return;
                catch
                    % fall through to pwelch
                end
            end

            % Pwelch fallback — matches computeASRMetrics.m exactly
            nseg = min(1024, max(256, 2^floor(log2(numel(x)))));
            win  = hamming(nseg);
            [Pxx, F] = pwelch(double(x), win, round(0.5*nseg), [], fs);
            mask = (F >= band(1)) & (F <= band(2));
            if any(mask)
                p = trapz(F(mask), Pxx(mask));
            else
                p = 0;
            end
        end

    end % static private methods

end
