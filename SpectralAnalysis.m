classdef SpectralAnalysis < handle
    % SpectralAnalysis - EEG band-power comparison before and after ASR cleaning.
    %
    %   OOP wrapper around the helperASR_metrics_* function chain.  All three
    %   helpers are faithfully replicated as static methods so no external
    %   helper files are required.
    %
    %   Helper → method mapping:
    %     helperASR_metrics_bandPower        → SpectralAnalysis.bandPower()
    %     helperASR_metrics_spectralAnalysis → called per segment inside compute()
    %     helperASR_metrics_spectralSummary  → results(k).summary.(segment)
    %
    %   ── PSD details ─────────────────────────────────────────────────────
    %   Welch PSD via pwelch(signal, window, noverlap, 0:freqResolution:fs/2, fs).
    %   Total power  = trapz(f, Pxx).
    %   Band power   = trapz over band-frequency indices.
    %   Relative power = bandPower / totalPower  (no eps guard — matches helper).
    %
    %   Configurable Welch parameters (set before compute(), or pass to compute()):
    %     welchWinSec   : window length in seconds   (default [] → MATLAB default)
    %     welchOverlap  : overlap fraction 0–<1       (default [] → MATLAB default)
    %     freqResolution: frequency axis step in Hz   (default 0.5 Hz)
    %
    %   Legacy default (welchWinSec=[], welchOverlap=[]) passes [] to pwelch for both
    %   window and noverlap, which is identical to the original helper behaviour.
    %
    %   ── Band definitions (from helperASR_metrics_spectralAnalysis) ──────
    %     delta  [1  4] Hz
    %     theta  [4  8] Hz
    %     alpha  [8 13] Hz
    %     beta  [13 30] Hz
    %     gamma [30 100] Hz
    %
    %   ── Result structure ────────────────────────────────────────────────
    %   results(k).perChannel.(seg).raw   — bandPower() output on raw signal
    %   results(k).perChannel.(seg).clean — bandPower() output on clean signal
    %   results(k).summary.(seg)          — channel-mean scalars (see below)
    %   results(k).psd.(seg).{raw_mean, clean_mean, f}  — mean PSD for plotting
    %   results(k).parameters             — parameter struct from subject_results
    %
    %   ── Summary struct layout ───────────────────────────────────────────
    %   Matches helperASR_metrics_spectralSummary output exactly:
    %     summary.(seg).<band>.totalPowerBefore
    %     summary.(seg).<band>.totalPowerAfter
    %     summary.(seg).<band>.relativePowerBefore
    %     summary.(seg).<band>.relativePowerAfter
    %     summary.(seg).totalPowerBefore
    %     summary.(seg).totalPowerAfter
    %
    %   Additional derived fields at the same level:
    %     summary.(seg).<band>.powerRatio     clean / raw absolute power
    %     summary.(seg).<band>.powerChange    (clean - raw) / raw
    %     summary.(seg).totalPowerChange      (clean - raw) / raw  for total
    %     summary.(seg).alphaOverAlphaBeta_raw    alpha/(alpha+beta)  raw
    %     summary.(seg).alphaOverAlphaBeta_clean  alpha/(alpha+beta)  clean
    %     summary.(seg).spectralCorr          Pearson corr of mean PSDs
    %     summary.(seg).spectralKL            KL(raw PSD || clean PSD)
    %
    %   ── getSummary key format ────────────────────────────────────────────
    %   Dot-delimited path into summary.(segment):
    %     getSummary('alpha.totalPowerAfter',   'closed')
    %     getSummary('beta.relativePowerBefore','open')
    %     getSummary('totalPowerChange',        'closed')
    %     getSummary('spectralCorr',            'closed')
    %
    %   ── Usage (via ExperimentAnalysis) ───────────────────────────────────
    %       ana.computeSpectral();
    %       ana.spectral.getSummary('alpha.totalPowerAfter', 'closed')
    %       ana.spectral.plotBands(1)
    %       ana.spectral.plotPSD(1, 'closed')
    %       T = ana.spectral.paramTable();
    %
    %   ── Usage (standalone) ───────────────────────────────────────────────
    %       sp = SpectralAnalysis(raw, subject_results, fs, subjectID, algorithmName);
    %       sp.compute();

    % ================================================================
    properties (Access = public)

        % Inputs
        raw               % struct: .calibration/.closed/.open [C x T]
        subject_results   % 1xN struct array (.cleanClosed/.cleanOpen/.cleanCalibration)
        fs                % Sampling rate Hz
        subjectID         % For display / titles only
        algorithmName     % For display / titles only
        nCombos

        % PSD / band config — may be overridden before compute(), or via compute() args
        freqResolution = 0.5   % Hz — frequency axis step for pwelch nfft argument
        welchWinSec    = []    % Welch window length (sec); [] = MATLAB default (~8× 1/freqRes)
        welchOverlap   = []    % Welch overlap fraction 0–<1; [] = MATLAB default (50%)

        bands = struct( ...    % Hz — matches helperASR_metrics_spectralAnalysis
            'delta', [1  4],  ...
            'theta', [4  8],  ...
            'alpha', [8  13], ...
            'beta',  [13 30], ...
            'gamma', [30 100])

        % Segments to analyse
        segments = {'closed', 'open', 'all'}

        % Optional streaming loader — set by ExperimentAnalysis for split-file layout.
        loaderFcn   = []   % @(k) ana.loadCombo(k)
        unloaderFcn = []   % @(k) ana.unloadCombo(k)

        % Results — 1xN struct array (populated by compute())
        results

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = SpectralAnalysis(raw, subject_results, fs, subjectID, algorithmName)
            % SPECTRALANALYSIS
            %   raw            : struct with .calibration/.closed/.open  [C x T]
            %   subject_results: 1xN struct array
            %                    (each needs .cleanCalibration/.cleanClosed/.cleanOpen)
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
            % COMPUTE - Run spectral analysis across all combos and segments.
            %
            %   Optional name-value options (override properties for this run):
            %     'welchWinSec'    scalar > 0   Welch window length (sec)
            %                                   [] = MATLAB default  (default [])
            %     'welchOverlap'   0 ≤ x < 1    Overlap fraction
            %                                   [] = MATLAB default 50%  (default [])
            %     'freqResolution' scalar > 0   Frequency axis step Hz  (default 0.5)
            %
            %   Example:
            %     sp.compute();                               % legacy defaults
            %     sp.compute('welchWinSec', 2, 'welchOverlap', 0.5);
            %     sp.compute('welchWinSec', 4);               % 4-sec Hann window

            p = inputParser;
            addParameter(p, 'welchWinSec',    obj.welchWinSec,    ...
                @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
            addParameter(p, 'welchOverlap',   obj.welchOverlap,   ...
                @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0 && x < 1));
            addParameter(p, 'freqResolution', obj.freqResolution, ...
                @(x) isnumeric(x) && isscalar(x) && x > 0);
            parse(p, varargin{:});

            obj.welchWinSec    = p.Results.welchWinSec;
            obj.welchOverlap   = p.Results.welchOverlap;
            obj.freqResolution = p.Results.freqResolution;

            % Resolve window / noverlap in samples once (shared by all pwelch calls)
            if isempty(obj.welchWinSec)
                welchWin = [];      % MATLAB default
            else
                welchWin = round(obj.welchWinSec * obj.fs);
            end
            if isempty(obj.welchOverlap) || isempty(welchWin)
                welchNov = [];      % MATLAB default (also safe if no explicit window)
            else
                welchNov = round(obj.welchOverlap * welchWin);
            end

            fprintf('SpectralAnalysis: S%d %s (%d combo(s))...\n', ...
                obj.subjectID, obj.algorithmName, obj.nCombos);

            % Raw segment signals — identical across all combos, build once
            rawMap = struct( ...
                'closed', double(obj.raw.closed), ...
                'open',   double(obj.raw.open),   ...
                'all',    double([obj.raw.calibration, obj.raw.closed, obj.raw.open]));

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

                cleanMap = struct( ...
                    'closed', double(sr.cleanClosed), ...
                    'open',   double(sr.cleanOpen),   ...
                    'all',    double([sr.cleanCalibration, sr.cleanClosed, sr.cleanOpen]));

                perCh  = struct();
                summ   = struct();
                psdOut = struct();

                for si = 1:numel(obj.segments)
                    seg = obj.segments{si};

                    Xr = rawMap.(seg);
                    Xc = cleanMap.(seg);

                    % Truncate to equal length (guard against carry off-by-one)
                    T  = min(size(Xr, 2), size(Xc, 2));
                    Xr = Xr(:, 1:T);
                    Xc = Xc(:, 1:T);

                    % helperASR_metrics_spectralAnalysis equivalent:
                    %   .before = bandPower(original)   .after = bandPower(cleaned)
                    raw_m   = SpectralAnalysis.bandPower(Xr, obj.fs, obj.bands, obj.freqResolution, welchWin, welchNov);
                    clean_m = SpectralAnalysis.bandPower(Xc, obj.fs, obj.bands, obj.freqResolution, welchWin, welchNov);

                    perCh.(seg).raw   = raw_m;
                    perCh.(seg).clean = clean_m;

                    % Mean PSD — needed for spectralCorr / spectralKL and plots
                    [psd_raw_mean,   f] = SpectralAnalysis.meanPSD(Xr, obj.fs, obj.freqResolution, welchWin, welchNov);
                    [psd_clean_mean, ~] = SpectralAnalysis.meanPSD(Xc, obj.fs, obj.freqResolution, welchWin, welchNov);

                    psdOut.(seg).raw_mean   = psd_raw_mean;
                    psdOut.(seg).clean_mean = psd_clean_mean;
                    psdOut.(seg).f          = f;

                    % helperASR_metrics_spectralSummary equivalent + derived extras
                    summ.(seg) = obj.buildSummary(raw_m, clean_m, psd_raw_mean, psd_clean_mean);
                end

                obj.results(c).perChannel = perCh;
                obj.results(c).summary    = summ;
                obj.results(c).psd        = psdOut;
                obj.results(c).parameters = sr.parameters;

                fprintf('  combo %d/%d done.\n', c, obj.nCombos);

                % Free combo signals if a loader is managing memory
                if ~isempty(obj.unloaderFcn)
                    obj.unloaderFcn(c);
                end
            end

            fprintf('  Done.\n');
        end

        % ── getSummary ───────────────────────────────────────────────
        function vals = getSummary(obj, fieldName, segment)
            % GETSUMMARY - Extract one scalar across all parameter combos.
            %
            %   fieldName : dot-delimited path into summary.(segment) struct.
            %               Band-level:   'alpha.totalPowerAfter'
            %               Top-level:    'totalPowerChange'
            %               Derived:      'spectralCorr'  |  'alphaOverAlphaBeta_clean'
            %   segment   : 'closed' (default) | 'open' | 'all'
            %
            %   Returns [1 x nCombos] double — NaN for missing entries.
            %
            %   Examples:
            %     sp.getSummary('alpha.totalPowerAfter',    'closed')
            %     sp.getSummary('delta.relativePowerBefore','open')
            %     sp.getSummary('gamma.powerRatio',         'closed')
            %     sp.getSummary('totalPowerChange',         'closed')
            %     sp.getSummary('spectralCorr',             'closed')
            %     sp.getSummary('alphaOverAlphaBeta_clean', 'closed')

            obj.checkComputed();
            if nargin < 3, segment = 'closed'; end

            vals  = NaN(1, obj.nCombos);
            parts = strsplit(fieldName, '.');

            for c = 1:obj.nCombos
                if ~isfield(obj.results(c).summary, segment), continue; end
                s = obj.results(c).summary.(segment);
                try
                    if numel(parts) == 1
                        vals(c) = s.(parts{1});
                    elseif numel(parts) == 2
                        vals(c) = s.(parts{1}).(parts{2});
                    end
                catch
                    % field absent for this combo — leave NaN
                end
            end
        end

        % ── paramTable ───────────────────────────────────────────────
        function T = paramTable(obj, segment)
            % PARAMTABLE - One row per combo: parameters + key spectral scalars.
            %   Band-level summary fields are flattened to column names like
            %   'alpha_pwBefore', 'alpha_relAfter', 'alpha_pwRatio', etc.
            obj.checkComputed();
            if nargin < 2, segment = 'closed'; end

            rows   = cell(obj.nCombos, 1);
            bnames = fieldnames(obj.bands);

            for c = 1:obj.nCombos
                row = table();
                row.comboIdx = c;

                % Parameter columns
                params = obj.results(c).parameters;
                pf = fieldnames(params);
                for i = 1:numel(pf)
                    val = params.(pf{i});
                    if isscalar(val) && isnumeric(val)
                        row.(pf{i}) = val;
                    end
                end

                % Spectral summary columns
                if isfield(obj.results(c).summary, segment)
                    s = obj.results(c).summary.(segment);

                    row.totalPowerBefore = s.totalPowerBefore;
                    row.totalPowerAfter  = s.totalPowerAfter;
                    row.totalPowerChange = s.totalPowerChange;
                    row.spectralCorr     = s.spectralCorr;
                    row.spectralKL       = s.spectralKL;
                    row.alphaOAB_raw     = s.alphaOverAlphaBeta_raw;
                    row.alphaOAB_clean   = s.alphaOverAlphaBeta_clean;

                    for i = 1:numel(bnames)
                        b = bnames{i};
                        row.(sprintf('%s_pwBefore',  b)) = s.(b).totalPowerBefore;
                        row.(sprintf('%s_pwAfter',   b)) = s.(b).totalPowerAfter;
                        row.(sprintf('%s_relBefore', b)) = s.(b).relativePowerBefore;
                        row.(sprintf('%s_relAfter',  b)) = s.(b).relativePowerAfter;
                        row.(sprintf('%s_pwRatio',   b)) = s.(b).powerRatio;
                        row.(sprintf('%s_pwChange',  b)) = s.(b).powerChange;
                    end
                end

                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

        % ── plotBands ────────────────────────────────────────────────
        function plotBands(obj, comboIdx, segment)
            % PLOTBANDS - Grouped bar: raw vs clean absolute band power.
            obj.checkComputed();
            if nargin < 2, comboIdx = 1; end
            if nargin < 3, segment  = 'closed'; end

            s      = obj.results(comboIdx).summary.(segment);
            bnames = fieldnames(obj.bands);
            nB     = numel(bnames);

            rawPow   = cellfun(@(b) s.(b).totalPowerBefore, bnames);
            cleanPow = cellfun(@(b) s.(b).totalPowerAfter,  bnames);

            figure('Name', sprintf('%s — band power S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));

            bh = bar(1:nB, [rawPow(:), cleanPow(:)], 'grouped');
            bh(1).FaceColor = [0.35 0.55 0.85];
            bh(2).FaceColor = [0.88 0.45 0.25];
            set(gca, 'XTick', 1:nB, 'XTickLabel', bnames);
            legend({'Raw', 'Clean'}, 'Location', 'northeast');
            ylabel('Mean band power (\muV^2/Hz \cdot s)');
            title(sprintf('%s — band power before/after ASR\nS%d  combo %d  (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter', 'none');
            grid on;
        end

        % ── plotBandChange ───────────────────────────────────────────
        function plotBandChange(obj, segment)
            % PLOTBANDCHANGE - Grouped bar of (clean-raw)/raw per band across combos.
            obj.checkComputed();
            if nargin < 2, segment = 'closed'; end

            bnames = fieldnames(obj.bands);
            nB     = numel(bnames);
            data   = zeros(obj.nCombos, nB);

            for c = 1:obj.nCombos
                s = obj.results(c).summary.(segment);
                for i = 1:nB
                    data(c, i) = s.(bnames{i}).powerChange;
                end
            end

            figure('Name', sprintf('%s — band power change (%s)', obj.algorithmName, segment));
            bar(data, 'grouped');
            set(gca, 'XTick', 1:obj.nCombos);
            xlabel('Parameter combo');
            ylabel('(Clean \minus Raw) / Raw');
            legend(bnames, 'Location', 'northeastoutside');
            title(sprintf('%s — fractional band power change (%s)', ...
                obj.algorithmName, segment), 'Interpreter', 'none');
            yline(0, 'k--', 'LineWidth', 1.2);
            grid on;
        end

        % ── plotPSD ──────────────────────────────────────────────────
        function plotPSD(obj, comboIdx, segment)
            % PLOTPSD - Log-scale mean PSD overlay: raw (blue) vs clean (red dashed).
            obj.checkComputed();
            if nargin < 2, comboIdx = 1; end
            if nargin < 3, segment  = 'closed'; end

            psd = obj.results(comboIdx).psd.(segment);

            figure('Name', sprintf('%s — PSD S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));

            semilogy(psd.f, psd.raw_mean,   'b-',  'LineWidth', 1.5, 'DisplayName', 'Raw');
            hold on;
            semilogy(psd.f, psd.clean_mean, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Clean');

            % Band boundary lines + labels
            bnames = fieldnames(obj.bands);
            edges  = unique([cellfun(@(b) obj.bands.(b)(1), bnames); ...
                             cellfun(@(b) obj.bands.(b)(2), bnames)]);
            for xi = 1:numel(edges)
                xline(edges(xi), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end

            ax = gca;  yl = ylim(ax);
            for i = 1:numel(bnames)
                xmid = mean(obj.bands.(bnames{i}));
                text(xmid, yl(2) * 0.65, bnames{i}, ...
                    'HorizontalAlignment', 'center', 'FontSize', 7, ...
                    'Color', [0.35 0.35 0.35], 'Interpreter', 'none');
            end

            legend('Location', 'northeast');
            xlabel('Frequency (Hz)');
            ylabel('Power (\muV^2/Hz)');
            xlim([0, min(100, obj.fs/2)]);
            title(sprintf('%s — mean PSD before/after ASR\nS%d  combo %d  (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter', 'none');
            grid on;
        end

        % ── plotRelativePower ────────────────────────────────────────
        function plotRelativePower(obj, comboIdx, segment)
            % PLOTRELATIVEPOWER - Stacked bar: relative band power before/after.
            obj.checkComputed();
            if nargin < 2, comboIdx = 1; end
            if nargin < 3, segment  = 'closed'; end

            s      = obj.results(comboIdx).summary.(segment);
            bnames = fieldnames(obj.bands);
            nB     = numel(bnames);

            rawRel   = cellfun(@(b) s.(b).relativePowerBefore, bnames);
            cleanRel = cellfun(@(b) s.(b).relativePowerAfter,  bnames);

            figure('Name', sprintf('%s — relative power S%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment));
            bar([rawRel(:)'; cleanRel(:)'], 'stacked');
            set(gca, 'XTick', 1:2, 'XTickLabel', {'Raw', 'Clean'});
            legend(bnames, 'Location', 'northeastoutside');
            ylabel('Relative power (fraction of total)');
            ylim([0 1]);  grid on;
            title(sprintf('%s — relative band power before/after ASR\nS%d combo %d (%s)', ...
                obj.algorithmName, obj.subjectID, comboIdx, segment), ...
                'Interpreter', 'none');
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function s = buildSummary(obj, raw_m, clean_m, psd_raw, psd_clean)
            % BUILDSUMMARY - Constructs the per-segment summary struct.
            %
            %   Top-level layout mirrors helperASR_metrics_spectralSummary exactly.
            %   Derived fields appended at the same level.

            bnames = fieldnames(obj.bands);

            % ── Exact replication of helperASR_metrics_spectralSummary ──
            for i = 1:numel(bnames)
                b = bnames{i};
                s.(b).totalPowerBefore    = mean(raw_m.bandPower.(b));
                s.(b).totalPowerAfter     = mean(clean_m.bandPower.(b));
                s.(b).relativePowerBefore = mean(raw_m.relativePower.(b));
                s.(b).relativePowerAfter  = mean(clean_m.relativePower.(b));
            end
            s.totalPowerBefore = mean(raw_m.total);
            s.totalPowerAfter  = mean(clean_m.total);

            % ── Derived: per-band power ratio and fractional change ──
            for i = 1:numel(bnames)
                b    = bnames{i};
                rPow = s.(b).totalPowerBefore;
                cPow = s.(b).totalPowerAfter;
                if rPow > 0
                    s.(b).powerRatio  = cPow / rPow;
                    s.(b).powerChange = (cPow - rPow) / rPow;
                else
                    s.(b).powerRatio  = NaN;
                    s.(b).powerChange = NaN;
                end
            end

            % ── Derived: total power fractional change ───────────────
            if s.totalPowerBefore > 0
                s.totalPowerChange = (s.totalPowerAfter - s.totalPowerBefore) / s.totalPowerBefore;
            else
                s.totalPowerChange = NaN;
            end

            % ── Derived: alpha concentration index ──────────────────
            % alpha / (alpha + beta) — sensitive to blink-band artefact removal
            s.alphaOverAlphaBeta_raw   = SpectralAnalysis.safeRatio( ...
                s.alpha.totalPowerBefore, s.alpha.totalPowerBefore + s.beta.totalPowerBefore);
            s.alphaOverAlphaBeta_clean = SpectralAnalysis.safeRatio( ...
                s.alpha.totalPowerAfter,  s.alpha.totalPowerAfter  + s.beta.totalPowerAfter);

            % ── Derived: spectral shape metrics ─────────────────────
            if ~isempty(psd_raw) && ~isempty(psd_clean)
                n  = min(numel(psd_raw), numel(psd_clean));
                pr = psd_raw(1:n);
                pc = psd_clean(1:n);

                % Pearson correlation of mean PSDs (spectral shape preservation)
                if std(pr) > 0 && std(pc) > 0
                    s.spectralCorr = corr(pr(:), pc(:));
                else
                    s.spectralCorr = NaN;
                end

                % KL divergence KL(raw || clean) normalised over frequency
                P = pr(:) / (sum(pr) + eps);
                Q = pc(:) / (sum(pc) + eps);
                P(P < 1e-12) = 1e-12;
                Q(Q < 1e-12) = 1e-12;
                s.spectralKL = sum(P .* log(P ./ Q));
            else
                s.spectralCorr = NaN;
                s.spectralKL   = NaN;
            end
        end

        function checkComputed(obj)
            if isempty(obj.results) || ~isfield(obj.results(1), 'summary')
                error('SpectralAnalysis:notComputed', ...
                    'Call compute() before accessing results.');
            end
        end

    end % private methods

    % ================================================================
    methods (Static, Access = private)

        function metrics = bandPower(eeg, fs, bands, freqRes, welchWin, welchNov)
            % BANDPOWER - Faithful replica of helperASR_metrics_bandPower.
            %
            %   Exact match to helper implementation:
            %     pwelch(eeg(ch,:), welchWin, welchNov, 0:fRes:fs/2, fs)
            %     relativePower = bandPow / totalPower   — no eps guard
            %     metrics.bandPower.(band)(ch,1) column assignment
            %
            %   Args:
            %     eeg       : [C x T] double
            %     fs        : sampling rate Hz
            %     bands     : struct of [lo hi] Hz band definitions
            %     freqRes   : frequency axis step Hz (default 0.5)
            %     welchWin  : pwelch window in samples  ([] = MATLAB default)
            %     welchNov  : pwelch noverlap in samples ([] = MATLAB default)

            if nargin < 4, freqRes  = 0.5; end
            if nargin < 5, welchWin = [];  end
            if nargin < 6, welchNov = [];  end

            [nCh, ~]  = size(eeg);
            bandNames = fieldnames(bands);

            metrics.total = zeros(nCh, 1);
            for i = 1:numel(bandNames)
                metrics.bandPower.(bandNames{i})     = zeros(nCh, 1);
                metrics.relativePower.(bandNames{i}) = zeros(nCh, 1);
            end

            for ch = 1:nCh
                [Pxx, f] = pwelch(double(eeg(ch, :)), welchWin, welchNov, 0:freqRes:fs/2, fs);

                totalPower        = trapz(f, Pxx);
                metrics.total(ch) = totalPower;

                for i = 1:numel(bandNames)
                    b    = bandNames{i};
                    idx  = f >= bands.(b)(1) & f < bands.(b)(2);
                    bPow = trapz(f(idx), Pxx(idx));

                    metrics.bandPower.(b)(ch, 1)     = bPow;
                    metrics.relativePower.(b)(ch, 1) = bPow / totalPower;
                end
            end
        end

        function [psd_mean, f] = meanPSD(eeg, fs, freqRes, welchWin, welchNov)
            % MEANPSD - Channel-mean Welch PSD (for spectralCorr/KL and plotting).
            %   Uses identical pwelch call to bandPower.
            %
            %   Args:
            %     eeg      : [C x T] double
            %     fs       : sampling rate Hz
            %     freqRes  : frequency axis step Hz (default 0.5)
            %     welchWin : pwelch window in samples  ([] = MATLAB default)
            %     welchNov : pwelch noverlap in samples ([] = MATLAB default)

            if nargin < 3, freqRes  = 0.5; end
            if nargin < 4, welchWin = [];  end
            if nargin < 5, welchNov = [];  end

            [nCh, ~] = size(eeg);
            psd_sum  = [];
            f        = [];

            for ch = 1:nCh
                [Pxx, fOut] = pwelch(double(eeg(ch, :)), welchWin, welchNov, 0:freqRes:fs/2, fs);
                if isempty(psd_sum)
                    psd_sum = zeros(size(Pxx));
                    f       = fOut;
                end
                psd_sum = psd_sum + Pxx;
            end

            psd_mean = psd_sum / nCh;
        end

        function r = safeRatio(num, denom)
            % SAFERATIO - num/denom, NaN if denom <= 0.
            if denom > 0
                r = num / denom;
            else
                r = NaN;
            end
        end

    end % static private methods

end
