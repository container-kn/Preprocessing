classdef probeAnalysis < handle
    % probeAnalysis - Unified post-hoc analyser for all ASR probe objects.
    %
    %   Accepts any asrProbe subclass (vanillaASRprobe, emaASRprobe,
    %   graphASRprobe, hmoASRprobe, hebbASRprobe) and computes a consistent
    %   summary of all time-series metrics the probe recorded.
    %
    %   Works on a 1xN array of subject_results structs (one per parameter
    %   combo) so it mirrors the BlinkAnalysis / SignalStatistics interface.
    %   Can also be constructed directly from a single probe object.
    %
    %   Usage (via ExperimentAnalysis — typical):
    %       ana = ExperimentAnalysis(folder, 'hmo-asr', 3, 500);
    %       ana.load();
    %       ana.computeProbe();       % creates ana.probeAna
    %       ana.probeAna.getSummary('trivialRate')
    %       ana.probeAna.getSummary('meanNormR')
    %       ana.probeAna.plotTimeSeries('normR', 1)   % combo 1
    %       ana.probeAna.summaryTable()
    %
    %   Usage (standalone, single probe):
    %       pa = probeAnalysis.fromProbe(myProbe, 'raw');
    %       pa.compute();
    %       T = pa.summaryTable();
    %
    %   Supported probe types (auto-detected from class name):
    %       asrProbe          — base metrics only
    %       vanillaASRprobe   — + subspaceAngle*, thresholdDrift
    %       emaASRprobe       — + subspaceAngle*, riemannDriftStep, thresholdDrift,
    %                             adaptationDeltaM/T, cleanBufferFill
    %       graphASRprobe     — + participationRatio, spectralCentroid,
    %                             lowHighEnergyRatio, riemannDriftStep
    %       hmoASRprobe       — + modelRiemannDrift, thresholdAdaptation,
    %                             sampleWeightTotal, basisCoherence,
    %                             riemannDistRaw, rejectionCount,
    %                             maxInflation, meanInflation, riemannDriftStep
    %       hebbASRprobe      — all of hmoASRprobe +
    %                             weightDeltaW, weightDeltaM, weightCovAlign,
    %                             mCondition, hebbStepCount
    %
    %   Summary scalars computed per time-series:
    %       mean, std, median, min, max, pct5, pct95, nNaN, nSamples
    %
    %   Additionally derived scalars (always present):
    %       trivialRate           — fraction of trivialFlag == 1
    %       meanRejectionCount    — mean(rejectionCount)  [hmo/hebb only]
    %       finalModelDrift       — last value of modelRiemannDrift  [hmo/hebb]
    %       hebbConvergenceStep   — first step where weightDeltaW < 0.01  [hebb]

    % ================================================================
    properties (Access = public)

        % Identity
        subjectID       = 0
        algorithmName   = 'unknown'
        probeRole       = 'raw'       % 'raw' | 'clean' — which probe stream
        nCombos         = 0

        % Input — either from subject_results or standalone
        subject_results               % 1xN struct array (may be empty in standalone mode)

        % Results — 1xN struct array
        %   results(k).scalars    : struct of summary scalars (all metrics)
        %   results(k).series     : struct of raw time-series (all probe fields)
        %   results(k).probeClass : string — detected probe class
        %   results(k).parameters : parameter struct from subject_results
        results

    end

    % ================================================================
    methods

        % ── Constructor (via subject_results — typical) ───────────────
        function obj = probeAnalysis(subject_results, probeRole, subjectID, algorithmName)
            % PROBEANALYSIS
            %   subject_results : 1xN struct array (each with .probeRaw/.probeClean)
            %   probeRole       : 'raw' | 'clean'  (default: 'raw')
            %   subjectID       : scalar (default: 0)
            %   algorithmName   : string (default: 'unknown')

            if nargin < 2 || isempty(probeRole),    probeRole    = 'raw';     end
            if nargin < 3 || isempty(subjectID),    subjectID    = 0;         end
            if nargin < 4 || isempty(algorithmName),algorithmName = 'unknown'; end

            if ~ismember(probeRole, {'raw','clean'})
                error('probeAnalysis:badRole', 'probeRole must be ''raw'' or ''clean''.');
            end

            obj.subject_results = subject_results;
            obj.probeRole       = probeRole;
            obj.subjectID       = subjectID;
            obj.algorithmName   = algorithmName;
            obj.nCombos         = numel(subject_results);
            obj.results         = struct();
        end

        % ── Static factory — standalone from a single probe ───────────
        function obj = fromProbe(probe, probeRole, algorithmName)
            % FROMPROBE - Build a probeAnalysis from a single probe object.
            if nargin < 2, probeRole    = 'raw';     end
            if nargin < 3, algorithmName = 'unknown'; end

            % Wrap the probe in a minimal subject_results struct
            sr = struct();
            if strcmp(probeRole, 'raw')
                sr.probeRaw   = probe;
                sr.probeClean = [];
            else
                sr.probeRaw   = [];
                sr.probeClean = probe;
            end
            sr.parameters = struct();

            obj = probeAnalysis(sr, probeRole, 0, algorithmName);
        end

        % ── Compute ──────────────────────────────────────────────────
        function compute(obj)
            % COMPUTE - Summarise probe time-series for all parameter combos.

            fprintf('probeAnalysis [%s]: S%d %s (%d combo(s))...\n', ...
                obj.probeRole, obj.subjectID, obj.algorithmName, obj.nCombos);

            for c = 1:obj.nCombos
                sr    = obj.subject_results(c);
                probe = obj.extractProbe(sr);

                if isempty(probe)
                    obj.results(c).scalars    = obj.nanScalars();
                    obj.results(c).series     = struct();
                    obj.results(c).probeClass = 'none';
                    if isfield(sr, 'parameters')
                        obj.results(c).parameters = sr.parameters;
                    else
                        obj.results(c).parameters = struct();
                    end
                    continue;
                end

                probeClass = class(probe);
                series     = obj.extractSeries(probe, probeClass);
                scalars    = obj.summariseSeries(series, probe, probeClass);

                obj.results(c).scalars    = scalars;
                obj.results(c).series     = series;
                obj.results(c).probeClass = probeClass;
                if isfield(sr, 'parameters')
                    obj.results(c).parameters = sr.parameters;
                else
                    obj.results(c).parameters = struct();
                end
            end

            fprintf('  Done.\n');
        end

        % ── getSummary ───────────────────────────────────────────────
        function vals = getSummary(obj, fieldName)
            % GETSUMMARY - Extract one scalar summary field across all combos.
            %   fieldName : any field in results(k).scalars
            %   Returns [1 x nCombos] vector (NaN for missing combos).
            %
            %   Common fields:
            %     'trivialRate'          'meanNormR'          'meanRiemannDrift'
            %     'meanCondCov'          'meanSpectralEntropy' 'meanEigInflationMax'
            %     'meanDirViolations'    'meanEnergyRatio'    'meanRejectionCount'
            %     'finalModelDrift'      'meanBasisCoherence' 'hebbConvergenceStep'
            %     'meanWeightDeltaW'     'meanWeightCovAlign' 'meanMCondition'

            obj.checkComputed();
            vals = NaN(1, obj.nCombos);
            for c = 1:obj.nCombos
                if isfield(obj.results(c).scalars, fieldName)
                    vals(c) = obj.results(c).scalars.(fieldName);
                end
            end
        end

        % ── summaryTable ─────────────────────────────────────────────
        function T = summaryTable(obj)
            % SUMMARYTABLE - One row per combo with parameters + all probe scalars.
            obj.checkComputed();

            rows = cell(obj.nCombos, 1);
            for c = 1:obj.nCombos
                row = table();
                row.comboIdx   = c;
                row.probeClass = string(obj.results(c).probeClass);

                % Parameters
                params = obj.results(c).parameters;
                f = fieldnames(params);
                for i = 1:numel(f)
                    val = params.(f{i});
                    if isscalar(val) && isnumeric(val)
                        row.(f{i}) = val;
                    end
                end

                % Scalar summary
                sc = obj.results(c).scalars;
                sf = fieldnames(sc);
                for i = 1:numel(sf)
                    val = sc.(sf{i});
                    if isscalar(val) && isnumeric(val)
                        row.(sf{i}) = val;
                    end
                end

                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

        % ── plotTimeSeries ───────────────────────────────────────────
        function plotTimeSeries(obj, fieldName, comboIdx)
            % PLOTTIMESERIES - Plot one probe time-series for a given combo.
            %   fieldName : name of a series field, e.g. 'normR', 'riemannDrift'
            %   comboIdx  : combo index (default: 1)
            obj.checkComputed();
            if nargin < 3, comboIdx = 1; end

            if comboIdx < 1 || comboIdx > obj.nCombos
                error('probeAnalysis:badCombo', 'comboIdx out of range.');
            end

            s = obj.results(comboIdx).series;
            if ~isfield(s, fieldName)
                fprintf('  [probeAnalysis] Field "%s" not present in combo %d.\n', ...
                    fieldName, comboIdx);
                return;
            end

            y = s.(fieldName);
            if isempty(y) || all(isnan(y))
                fprintf('  [probeAnalysis] Field "%s" is empty/all-NaN in combo %d.\n', ...
                    fieldName, comboIdx);
                return;
            end

            figure('Name', sprintf('%s: %s (combo %d)', obj.algorithmName, fieldName, comboIdx));
            plot(y, '-', 'LineWidth', 1.2);
            xlabel('Window index');
            ylabel(fieldName, 'Interpreter', 'none');
            title(sprintf('%s / %s — S%d combo %d', ...
                obj.algorithmName, fieldName, obj.subjectID, comboIdx), ...
                'Interpreter', 'none');
            grid on;
        end

        % ── plotAllSeries ────────────────────────────────────────────
        function plotAllSeries(obj, comboIdx)
            % PLOTALLSERIES - Subplot grid of all non-trivial series for one combo.
            obj.checkComputed();
            if nargin < 2, comboIdx = 1; end

            s = obj.results(comboIdx).series;
            fields = fieldnames(s);

            % Filter to non-empty numeric arrays
            valid = {};
            for i = 1:numel(fields)
                v = s.(fields{i});
                if isnumeric(v) && ~isempty(v) && ~all(isnan(v))
                    valid{end+1} = fields{i}; %#ok<AGROW>
                end
            end

            if isempty(valid)
                fprintf('  [probeAnalysis] No plottable series in combo %d.\n', comboIdx);
                return;
            end

            nc = ceil(sqrt(numel(valid)));
            nr = ceil(numel(valid) / nc);
            figure('Name', sprintf('%s probe series (combo %d)', obj.algorithmName, comboIdx));
            for i = 1:numel(valid)
                subplot(nr, nc, i);
                plot(s.(valid{i}), '-', 'LineWidth', 1);
                title(valid{i}, 'Interpreter', 'none', 'FontSize', 8);
                xlabel('win');  grid on;
            end
            sgtitle(sprintf('%s / %s — combo %d', ...
                obj.algorithmName, obj.probeRole, comboIdx), 'Interpreter', 'none');
        end

    end % public methods

    % ================================================================
    methods (Static)

        % ── Static factory — accessible without instance ──────────────
        function obj = fromProbeStatic(probe, probeRole, algorithmName)
            % FROMPROBSTATIC - Build probeAnalysis from a single probe.
            if nargin < 2, probeRole    = 'raw';     end
            if nargin < 3, algorithmName = 'unknown'; end
            sr = struct('parameters', struct());
            if strcmp(probeRole, 'raw')
                sr.probeRaw = probe;  sr.probeClean = [];
            else
                sr.probeRaw = [];  sr.probeClean = probe;
            end
            obj = probeAnalysis(sr, probeRole, 0, algorithmName);
        end

    end % static methods

    % ================================================================
    methods (Access = private)

        % ── extractProbe ─────────────────────────────────────────────
        function probe = extractProbe(obj, sr)
            % EXTRACTPROBE - Pull the appropriate probe from a subject_results entry.
            probe = [];
            field = sprintf('probe%s', [upper(obj.probeRole(1)), obj.probeRole(2:end)]);
            if isfield(sr, field) && ~isempty(sr.(field))
                probe = sr.(field);
                return;
            end
            % emaASR stores its probe in .emaProbe
            if strcmp(obj.probeRole, 'raw') && isfield(sr, 'emaProbe') && ~isempty(sr.emaProbe)
                probe = sr.emaProbe;
            end
        end

        % ── extractSeries ────────────────────────────────────────────
        function series = extractSeries(~, probe, probeClass)
            % EXTRACTSERIES - Pull all numeric time-series from a probe.
            %   Returns a struct with one field per series.

            % ---- Fields shared by ALL probe types (asrProbe base) ----
            baseFields = { ...
                'traceCov', 'maxEigCov', 'condCov', ...
                'spectralEntropy', 'logDetDrift', 'riemannDrift', ...
                'eigInflationMax', 'eigInflationTopK', ...
                'dirInflationMax', 'dirInflationMean', 'dirNumViolations', ...
                'energyRatio', 'normR', 'rankR', 'trivialFlag' };

            % ---- Fields by probe subclass ----
            extraFields = containers.Map();
            extraFields('vanillaASRprobe') = { ...
                'subspaceAngleMax', 'subspaceAngleMean', 'thresholdDrift' };
            extraFields('emaASRprobe') = { ...
                'riemannDriftStep', 'thresholdDrift', ...
                'subspaceAngleMax', 'subspaceAngleMean', 'subspaceAngleStep', ...
                'adaptationDeltaM', 'adaptationDeltaT', 'cleanBufferFill' };
            extraFields('graphASRprobe') = { ...
                'participationRatio', 'spectralCentroid', 'lowHighEnergyRatio', ...
                'riemannDriftStep', 'spectralDisagreement', 'subspaceAngle', 'commutatorNorm' };
            extraFields('hmoASRprobe') = { ...
                'modelRiemannDrift', 'thresholdAdaptation', 'sampleWeightTotal', ...
                'basisCoherence', 'riemannDistRaw', 'rejectionCount', ...
                'maxInflation', 'meanInflation', 'riemannDriftStep' };
            extraFields('hebbASRprobe') = { ...
                'modelRiemannDrift', 'thresholdAdaptation', 'sampleWeightTotal', ...
                'basisCoherence', 'riemannDistRaw', 'rejectionCount', ...
                'maxInflation', 'meanInflation', 'riemannDriftStep', ...
                'weightDeltaW', 'weightDeltaM', 'weightCovAlign', ...
                'mCondition', 'hebbStepCount' };

            % Collect all relevant field names for this probe class
            allFields = baseFields;
            % Walk up the inheritance: hebbASR has all hmo fields too
            if isKey(extraFields, probeClass)
                allFields = [allFields, extraFields(probeClass)];
            end

            series = struct();
            for i = 1:numel(allFields)
                fn = allFields{i};
                if isprop(probe, fn) && isnumeric(probe.(fn)) && ~isempty(probe.(fn))
                    series.(fn) = probe.(fn)(:);   % always column vector
                end
            end
        end

        % ── summariseSeries ──────────────────────────────────────────
        function sc = summariseSeries(~, series, probe, probeClass)
            % SUMMARISESERIES - Compute summary scalars from all series.

            sc = struct();
            fields = fieldnames(series);

            % --- Per-series summary stats ---
            for i = 1:numel(fields)
                fn = fields{i};
                v  = series.(fn);
                v  = double(v(isfinite(v)));    % strip NaN/Inf for stats

                prefix = sprintf('mean%s%s', upper(fn(1)), fn(2:end));

                if isempty(v)
                    sc.(sprintf('mean_%s', fn))   = NaN;
                    sc.(sprintf('std_%s',  fn))   = NaN;
                    sc.(sprintf('med_%s',  fn))   = NaN;
                    sc.(sprintf('max_%s',  fn))   = NaN;
                    sc.(sprintf('p95_%s',  fn))   = NaN;
                else
                    sc.(sprintf('mean_%s', fn))   = mean(v);
                    sc.(sprintf('std_%s',  fn))   = std(v);
                    sc.(sprintf('med_%s',  fn))   = median(v);
                    sc.(sprintf('max_%s',  fn))   = max(v);
                    sc.(sprintf('p95_%s',  fn))   = prctile(v, 95);
                end
                % Convenient camelCase alias for getSummary
                sc.(prefix) = sc.(sprintf('mean_%s', fn));
            end

            % --- Derived scalars ---

            % trivialRate
            if isfield(series, 'trivialFlag')
                sc.trivialRate = mean(series.trivialFlag == 1, 'omitnan');
            else
                sc.trivialRate = NaN;
            end

            % dirViolationRate
            if isfield(series, 'dirNumViolations')
                sc.meanDirViolations = sc.(sprintf('mean_%s', 'dirNumViolations'));
            end

            % hmoASR / hebbASR extras
            if isfield(series, 'rejectionCount')
                sc.meanRejectionCount = sc.mean_rejectionCount;
                sc.maxRejectionCount  = sc.max_rejectionCount;
            else
                sc.meanRejectionCount = NaN;
                sc.maxRejectionCount  = NaN;
            end

            if isfield(series, 'modelRiemannDrift') && ~isempty(series.modelRiemannDrift)
                v = series.modelRiemannDrift;
                sc.finalModelDrift = v(end);
                sc.maxModelDrift   = max(v);
            else
                sc.finalModelDrift = NaN;
                sc.maxModelDrift   = NaN;
            end

            if isfield(series, 'basisCoherence')
                sc.meanBasisCoherence = sc.mean_basisCoherence;
            else
                sc.meanBasisCoherence = NaN;
            end

            % hebbASR convergence step
            if isfield(series, 'weightDeltaW')
                dW  = series.weightDeltaW;
                idx = find(dW < 0.01, 1, 'first');
                sc.hebbConvergenceStep = idx;   % NaN if never converged
                if isempty(sc.hebbConvergenceStep)
                    sc.hebbConvergenceStep = NaN;
                end
                sc.meanWeightDeltaW   = sc.mean_weightDeltaW;
                sc.meanWeightCovAlign = sc.(sprintf('mean_%s', 'weightCovAlign'));
            else
                sc.hebbConvergenceStep = NaN;
                sc.meanWeightDeltaW    = NaN;
                sc.meanWeightCovAlign  = NaN;
            end

            if isfield(series, 'mCondition')
                sc.meanMCondition = sc.mean_mCondition;
                sc.maxMCondition  = sc.max_mCondition;
            else
                sc.meanMCondition = NaN;
                sc.maxMCondition  = NaN;
            end

            % Convenient aliases matching SignalStatistics.computeProbeStats names
            sc.meanNormR        = sc.(sprintf('mean_%s', 'normR'));
            sc.meanRiemannDrift = sc.(sprintf('mean_%s', 'riemannDrift'));
            sc.meanCondCov      = sc.(sprintf('mean_%s', 'condCov'));
            sc.meanSpectralEntropy = sc.(sprintf('mean_%s', 'spectralEntropy'));

            % Total window count
            sc.nWindows = probe.step;
        end

        % ── nanScalars ───────────────────────────────────────────────
        function sc = nanScalars(~)
            % NANSCALARS - Return a scalar struct of NaN for missing probes.
            sc = struct( ...
                'trivialRate',         NaN, ...
                'meanNormR',           NaN, ...
                'meanRiemannDrift',    NaN, ...
                'meanCondCov',         NaN, ...
                'meanSpectralEntropy', NaN, ...
                'meanRejectionCount',  NaN, ...
                'finalModelDrift',     NaN, ...
                'meanBasisCoherence',  NaN, ...
                'hebbConvergenceStep', NaN, ...
                'meanWeightDeltaW',    NaN, ...
                'meanWeightCovAlign',  NaN, ...
                'meanMCondition',      NaN, ...
                'nWindows',            0);
        end

        % ── checkComputed ────────────────────────────────────────────
        function checkComputed(obj)
            if isempty(obj.results) || ~isfield(obj.results(1), 'scalars')
                error('probeAnalysis:notComputed', 'Call compute() before accessing results.');
            end
        end

    end % private methods

end
