classdef SignalStatistics < handle
    % SignalStatistics - Descriptive and reconstruction statistics for ASR results.
    %
    %   Computes signal-level statistics across three groupings:
    %
    %   GROUP 1 — Cognitive state / segment
    %     .calibration, .closed, .open, .all
    %     Each covers the full raw and cleaned signal for that segment.
    %
    %   GROUP 2 — Reconstruction activity (from modifiedMask)
    %     .trivial    — samples where ASR passed through unchanged
    %     .nontrivial — samples where reconstruction was triggered
    %
    %   Statistics computed per group on raw and clean signals:
    %     variance, std, RMS, peak-to-peak amplitude, kurtosis, skewness
    %     (all stored per-channel; scalar summary = cross-channel mean)
    %
    %   Statistics computed on diff (clean - raw) signal:
    %     signalChange          ||clean-raw||_F / ||raw||_F      global scalar
    %     snrDB                 10*log10(power_raw/power_diff)   mean across channels
    %     rrmse                 mean per-channel RRMSE           mean(rms(diff)/rms(raw))
    %     corr                  mean per-channel Pearson corr    mean(corr(raw,clean))
    %     rrmse_per_channel     [C x 1]  per-channel RRMSE
    %     corr_per_channel      [C x 1]  per-channel Pearson corr
    %
    %   Note: rrmse and corr are whole-signal metrics.  For the same metrics
    %   restricted to blink-free regions, see QuietRegionAnalysis.
    %
    %   Additionally from the probe (per combo, if available):
    %     trivialRate      — fraction of windows that were trivial
    %     meanNormR        — mean reconstruction matrix distance from identity
    %     meanRiemannDrift — mean AIRM distance from calibration baseline
    %
    %   Usage (via ExperimentAnalysis — typical):
    %       ana.computeSignalStats();
    %       ana.signalStats.results(k).segment.closed.diff.rrmse
    %       ana.signalStats.getSummary('rrmse', 'closed', 'diff')
    %       ana.signalStats.getSummary('corr',  'closed', 'diff')
    %       ana.signalStats.getSummary('signalChange', 'closed', 'diff')
    %
    %   Usage (standalone):
    %       ss = SignalStatistics(raw, subject_results, fs, subjectID, algorithmName);
    %       ss.compute();

    properties (Access = public)
        % Inputs
        raw
        subject_results
        fs
        subjectID
        algorithmName
        nCombos

        % Results — 1xN struct array
        %   results(k).segment.(calibration|closed|open|all)
        %                 .raw   : stats on raw signal
        %                 .clean : stats on cleaned signal
        %                 .diff  : stats on (clean - raw) difference
        %   results(k).window.(trivial|nontrivial)
        %                 .raw   : stats on raw samples in those windows
        %                 .clean : stats on cleaned samples
        %   results(k).probe    : stats derived from probe telemetry
        results
    end

    methods

        % ================================================================
        % Constructor
        % ================================================================
        function obj = SignalStatistics(raw, subject_results, fs, subjectID, algorithmName)
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

        % ================================================================
        % compute
        % ================================================================
        function compute(obj)
            % COMPUTE - Run all statistics across all combos.

            fprintf('SignalStatistics: S%d %s (%d combo(s))...\n', ...
                obj.subjectID, obj.algorithmName, obj.nCombos);

            % Raw segments — identical across all combos, compute once
            rawSegs = struct( ...
                'calibration', obj.raw.calibration, ...
                'closed',      obj.raw.closed,      ...
                'open',        obj.raw.open,        ...
                'all',         [obj.raw.calibration, obj.raw.closed, obj.raw.open]);

            for c = 1:obj.nCombos
                sr = obj.subject_results(c);

                % ---- GROUP 1: per-segment stats ----
                cleanSegs = struct( ...
                    'calibration', sr.cleanCalibration, ...
                    'closed',      sr.cleanClosed,      ...
                    'open',        sr.cleanOpen,        ...
                    'all',         [sr.cleanCalibration, sr.cleanClosed, sr.cleanOpen]);

                segNames = {'calibration','closed','open','all'};
                segStats = struct();
                for s = 1:numel(segNames)
                    nm = segNames{s};
                    Xr = double(rawSegs.(nm));
                    Xc = double(cleanSegs.(nm));
                    segStats.(nm).raw   = SignalStatistics.computeStats(Xr);
                    segStats.(nm).clean = SignalStatistics.computeStats(Xc);
                    segStats.(nm).diff  = SignalStatistics.computeDiff(Xr, Xc);
                end

                % ---- GROUP 2: trivial / non-trivial window stats ----
                windowStats = obj.computeWindowStats(sr, rawSegs.all, cleanSegs.all);

                % ---- Probe-derived summary ----
                probeStats = obj.computeProbeStats(sr);

                obj.results(c).segment = segStats;
                obj.results(c).window  = windowStats;
                obj.results(c).probe   = probeStats;
            end

            fprintf('  Done.\n');
        end

        % ================================================================
        % getSummary — extract one scalar field across all combos
        % ================================================================
        function vals = getSummary(obj, fieldPath, segment, signalType)
            % GETSUMMARY - Pull one scalar across all combos.
            %
            %   fieldPath  : stat field name e.g. 'rms', 'kurtosis', 'trivialRate'
            %   segment    : 'calibration'|'closed'|'open'|'all'  (default: 'closed')
            %                or 'trivial'|'nontrivial' for window group
            %                or 'probe' for probe-derived fields
            %   signalType : 'raw'|'clean'|'diff'  (default: 'clean')
            %                (ignored for probe and window group)
            %
            %   Examples:
            %       ss.getSummary('rms',          'closed', 'clean')
            %       ss.getSummary('trivialRate',  'probe')
            %       ss.getSummary('signalChange', 'all',    'diff')

            obj.checkComputed();
            if nargin < 3, segment    = 'closed'; end
            if nargin < 4, signalType = 'clean';  end

            vals = zeros(1, obj.nCombos);
            for c = 1:obj.nCombos
                if strcmp(segment, 'probe')
                    vals(c) = obj.results(c).probe.(fieldPath);
                elseif ismember(segment, {'trivial','nontrivial'})
                    vals(c) = obj.results(c).window.(segment).(signalType).(fieldPath);
                else
                    vals(c) = obj.results(c).segment.(segment).(signalType).(fieldPath);
                end
            end
        end

        % ================================================================
        % paramTable
        % ================================================================
        function T = paramTable(obj)
            % PARAMTABLE - One row per combo with parameters + key summary stats.
            obj.checkComputed();
            rows = cell(obj.nCombos, 1);
            for c = 1:obj.nCombos
                params = obj.subject_results(c).parameters;
                fields = fieldnames(params);
                row    = table();
                row.comboIdx = c;
                for f = 1:numel(fields)
                    val = params.(fields{f});
                    if isscalar(val) && isnumeric(val)
                        row.(fields{f}) = val;
                    end
                end
                % Attach key summary stats
                row.trivialRate_closed  = obj.results(c).segment.closed.clean.trivialRate;
                row.rms_clean_closed    = obj.results(c).segment.closed.clean.rms;
                row.rms_raw_closed      = obj.results(c).segment.closed.raw.rms;
                row.signalChange_closed = obj.results(c).segment.closed.diff.signalChange;
                if isfield(obj.results(c).probe, 'trivialRate')
                    row.probe_trivialRate = obj.results(c).probe.trivialRate;
                end
                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

    end % public methods

    % ================================================================
    % Private helpers
    % ================================================================
    methods (Access = private)

        function ws = computeWindowStats(obj, sr, Xraw_all, Xclean_all)
            % COMPUTEWINDOWSTATS - Split samples by trivial/nontrivial using modifiedMask.
            %   modifiedMask is [1 x T_all] where T_all = sum of all three segments.
            %   true  = reconstruction happened (non-trivial)
            %   false = passed through unchanged (trivial)

            ws = struct();

            mask = sr.modifiedMask(:)';   % ensure row vector

            % Guard: mask might be shorter than signal if processing was partial
            T = min(numel(mask), size(Xraw_all, 2));
            if T == 0
                ws.trivial    = SignalStatistics.emptyWindowStats();
                ws.nontrivial = SignalStatistics.emptyWindowStats();
                return;
            end

            mask       = mask(1:T);
            Xr         = Xraw_all(:,   1:T);
            Xc         = Xclean_all(:, 1:T);

            trivIdx    = ~mask;
            nontrivIdx = mask;

            if any(trivIdx)
                ws.trivial.raw   = SignalStatistics.computeStats(Xr(:, trivIdx));
                ws.trivial.clean = SignalStatistics.computeStats(Xc(:, trivIdx));
                ws.trivial.diff  = SignalStatistics.computeDiff( Xr(:, trivIdx), Xc(:, trivIdx));
            else
                ws.trivial = SignalStatistics.emptyWindowStats();
            end

            if any(nontrivIdx)
                ws.nontrivial.raw   = SignalStatistics.computeStats(Xr(:, nontrivIdx));
                ws.nontrivial.clean = SignalStatistics.computeStats(Xc(:, nontrivIdx));
                ws.nontrivial.diff  = SignalStatistics.computeDiff( Xr(:, nontrivIdx), Xc(:, nontrivIdx));
            else
                ws.nontrivial = SignalStatistics.emptyWindowStats();
            end

            % Attach trivial rate as a convenience scalar to each group
            trivRate = mean(trivIdx);
            ws.trivial.trivialRate    = trivRate;
            ws.nontrivial.trivialRate = trivRate;   % same scalar, for convenience
        end

        function ps = computeProbeStats(~, sr)
            % COMPUTEPROBSTATS - Summary scalars from the probe object.
            %   Probes may be absent (e.g. sweep-mode accumulated files where
            %   probes weren't saved) — all fields fall back to NaN gracefully.

            ps = struct( ...
                'trivialRate',      NaN, ...
                'meanNormR',        NaN, ...
                'meanRiemannDrift', NaN, ...
                'meanCondCov',      NaN, ...
                'meanSpectralEntropy', NaN);

            probe = [];
            if isfield(sr, 'probeRaw') && ~isempty(sr.probeRaw)
                probe = sr.probeRaw;
            end
            if isempty(probe), return; end

            if isprop(probe, 'trivialFlag') && ~isempty(probe.trivialFlag)
                ps.trivialRate = mean(probe.trivialFlag, 'omitnan');
            end
            if isprop(probe, 'normR') && ~isempty(probe.normR)
                ps.meanNormR = mean(probe.normR, 'omitnan');
            end
            if isprop(probe, 'riemannDrift') && ~isempty(probe.riemannDrift)
                ps.meanRiemannDrift = mean(probe.riemannDrift, 'omitnan');
            end
            if isprop(probe, 'condCov') && ~isempty(probe.condCov)
                ps.meanCondCov = mean(probe.condCov, 'omitnan');
            end
            if isprop(probe, 'spectralEntropy') && ~isempty(probe.spectralEntropy)
                ps.meanSpectralEntropy = mean(probe.spectralEntropy, 'omitnan');
            end
        end

        function checkComputed(obj)
            if isempty(obj.results) || ~isfield(obj.results(1), 'segment')
                error('SignalStatistics:notComputed', 'Call compute() before accessing results.');
            end
        end

    end % private methods

    % ================================================================
    % Static helpers
    % ================================================================
    methods (Static, Access = private)

        function s = computeStats(X)
            % COMPUTESTATS - Per-channel + cross-channel mean statistics.
            %   X : [C x T] signal matrix

            if isempty(X)
                s = SignalStatistics.nanStats();
                return;
            end

            X = double(X);

            % Per-channel stats [C x 1]
            s.variance   = var(X, 0, 2);
            s.std        = std(X, 0, 2);
            s.rms        = sqrt(mean(X.^2, 2));
            s.peakToPeak = max(X,[],2) - min(X,[],2);
            s.kurtosis   = kurtosis(X, 1, 2);
            s.skewness   = skewness(X, 1, 2);
            s.mean       = mean(X, 2);

            % Cross-channel means (scalars — convenient for sweeps)
            s.rms_mean        = mean(s.rms);
            s.variance_mean   = mean(s.variance);
            s.kurtosis_mean   = mean(s.kurtosis);
            s.skewness_mean   = mean(s.skewness);
            s.peakToPeak_mean = mean(s.peakToPeak);

            % Kept as scalar aliases for getSummary convenience
            s.rms       = s.rms_mean;
            s.kurtosis  = s.kurtosis_mean;
            s.skewness  = s.skewness_mean;
            % (perChannel versions available via results directly)
        end

        function d = computeDiff(Xraw, Xclean)
            % COMPUTEDIFF - Statistics of (clean - raw) difference signal.
            %
            %   Scalar summary fields (cross-channel means, accessible via getSummary):
            %     signalChange   ||clean-raw||_F / ||raw||_F       global Frobenius ratio
            %     snrDB          10*log10(power_raw/power_diff)     mean across channels
            %     rrmse          mean per-channel RRMSE             mean(rms(diff)/rms(raw))
            %     corr           mean per-channel Pearson corr      mean(corr(raw,clean))
            %
            %   Per-channel arrays (accessible via results(k).segment directly):
            %     rrmse_per_channel   [C x 1]  rms(diff) / rms(raw)
            %     corr_per_channel    [C x 1]  Pearson corr(raw, clean)
            %
            %   Relationship to other classes:
            %     QuietRegionAnalysis.rre / .corr — same metrics on QUIET samples only.
            %     These fields are whole-signal equivalents for comparison.

            if isempty(Xraw) || isempty(Xclean)
                d = SignalStatistics.nanDiffStats();
                return;
            end

            Xraw   = double(Xraw);
            Xclean = double(Xclean);
            T      = min(size(Xraw,2), size(Xclean,2));
            Xraw   = Xraw(:,   1:T);
            Xclean = Xclean(:, 1:T);
            diff   = Xclean - Xraw;
            nCh    = size(Xraw, 1);

            d = SignalStatistics.computeStats(diff);

            % --- Global signal change: ||clean-raw||_F / ||raw||_F ---
            normRaw = norm(Xraw, 'fro');
            if normRaw > 0
                d.signalChange = norm(diff, 'fro') / normRaw;
            else
                d.signalChange = NaN;
            end

            % --- SNR in dB: 10*log10(power_raw / power_diff), mean across channels ---
            power_raw  = mean(Xraw.^2,  2);
            power_diff = mean(diff.^2,  2);
            snr_per_ch = 10 * log10(power_raw ./ max(power_diff, eps));
            d.snrDB    = mean(snr_per_ch, 'omitnan');

            % --- Per-channel RRMSE: rms(diff) / rms(raw) ---
            % Equivalent to helperASR_metrics_RRMSE but without the +eps denominator guard
            % (we use max(...,eps) to avoid division by zero on flat channels).
            rms_raw  = sqrt(mean(Xraw.^2,  2));   % [C x 1]
            rms_diff = sqrt(mean(diff.^2,  2));   % [C x 1]
            rrmse_ch = rms_diff ./ max(rms_raw, eps);
            d.rrmse_per_channel = rrmse_ch;
            d.rrmse             = mean(rrmse_ch, 'omitnan');

            % --- Per-channel Pearson correlation: corr(raw, clean) ---
            corr_ch = NaN(nCh, 1);
            for ch = 1:nCh
                xr = Xraw(ch, :);
                xc = Xclean(ch, :);
                if std(xr) > 0 && std(xc) > 0
                    c = corrcoef(xr, xc);
                    corr_ch(ch) = c(1,2);
                end
            end
            d.corr_per_channel = corr_ch;
            d.corr             = mean(corr_ch, 'omitnan');
        end

        function s = nanStats()
            s = struct( ...
                'variance', NaN, 'std', NaN, 'rms', NaN, ...
                'peakToPeak', NaN, 'kurtosis', NaN, 'skewness', NaN, ...
                'mean', NaN, 'rms_mean', NaN, 'variance_mean', NaN, ...
                'kurtosis_mean', NaN, 'skewness_mean', NaN, ...
                'peakToPeak_mean', NaN);
        end

        function s = nanDiffStats()
            % NANDIFFSTATS - Full NaN struct for computeDiff output (adds diff-specific fields).
            s = SignalStatistics.nanStats();
            s.signalChange      = NaN;
            s.snrDB             = NaN;
            s.rrmse             = NaN;
            s.corr              = NaN;
            s.rrmse_per_channel = NaN;
            s.corr_per_channel  = NaN;
        end

        function ws = emptyWindowStats()
            nan_s  = SignalStatistics.nanStats();
            nan_ds = SignalStatistics.nanDiffStats();
            ws.raw        = nan_s;
            ws.clean      = nan_s;
            ws.diff       = nan_ds;
            ws.trivialRate = NaN;
        end

    end % static methods

end
