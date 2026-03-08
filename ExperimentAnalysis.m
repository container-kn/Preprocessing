classdef ExperimentAnalysis < handle
    % ExperimentAnalysis - Data loader and analyser registry for ASR results.
    %
    %   Responsible for loading data (three modes) and instantiating metric
    %   analysers. Each metric lives in its own class — ExperimentAnalysis
    %   is purely a container and orchestrator.
    %
    %   THREE LOAD MODES:
    %
    %   SWEEP MODE (default constructor) — reads accumulated file:
    %       ana = ExperimentAnalysis(accumulatedFolder, 'vanilla-asr', 3, 500);
    %       ana.load();
    %
    %   SINGLE RUN MODE — reads one individual result .mat:
    %       ana = ExperimentAnalysis.fromSingleRun(resultsFolder, 'vanilla-asr', 3, 500, params);
    %       ana.load();
    %
    %   MEMORY MODE — no file I/O, data from experiment.outputResults():
    %       [data, asrObj] = experiment.outputResults();
    %       ana = ExperimentAnalysis.fromOutputResults(data, asrObj, 'vanilla-asr', 500);
    %
    %   RUNNING METRICS:
    %       ana.computeBlinks();              % delegates to BlinkAnalysis
    %       % future: ana.computeBandPower();
    %       % future: ana.computeSNR();
    %
    %   ACCESSING RESULTS:
    %       ana.blink                         % BlinkAnalysis instance
    %       ana.blink.results(k).summary
    %       ana.blink.results(k).perChannel
    %       ana.blink.getSummary('blinkReductionRatio_closed')
    %       ana.blink.paramTable()
    %
    %   BACKWARD-COMPATIBLE SHORTHAND (mirrors old metrics struct layout):
    %       ana.metrics(k).blink   → ana.blink.results(k)

    properties (Access = public)
        % Identity
        accumulatedFolder
        resultsFolder
        algorithmName
        subjectID
        fs
        mode              % 'sweep' | 'single' | 'memory'

        % Loaded data — shared across all analysers
        raw               % struct: .calibration/.closed/.open [C x T]
        subject_results   % 1xN struct array
        nCombos

        % ---- Analyser instances (populated by compute* methods) ----
        blink             % BlinkAnalysis instance
        signalStats       % SignalStatistics instance
        ica               % ICAAnalysis instance
        % bandPower       % BandPowerAnalysis instance  (future)
        % snr             % SNRAnalysis instance        (future)
    end

    properties (Access = private)
        singleParams
    end

    methods

        % ================================================================
        % Constructor — sweep mode
        % ================================================================
        function obj = ExperimentAnalysis(accumulatedFolder, algorithmName, subjectID, fs)
            if ~isfolder(accumulatedFolder)
                error('ExperimentAnalysis:invalidFolder', ...
                    'accumulatedFolder does not exist: %s', accumulatedFolder);
            end
            obj.accumulatedFolder = accumulatedFolder;
            obj.algorithmName     = algorithmName;
            obj.subjectID         = subjectID;
            obj.fs                = fs;
            obj.mode              = 'sweep';
        end

        % ================================================================
        % Load — dispatches on mode
        % ================================================================
        function load(obj)
            switch obj.mode
                case 'sweep',  obj.loadAccumulated();
                case 'single', obj.loadSingleRun();
                case 'memory', fprintf('Memory mode -- data already loaded.\n');
                otherwise
                    error('ExperimentAnalysis:unknownMode', 'Unknown mode: %s', obj.mode);
            end
        end

        % ================================================================
        % computeBlinks — delegates to BlinkAnalysis
        % ================================================================
        function computeBlinks(obj, varargin)
            % COMPUTEBLINKS - Instantiate BlinkAnalysis and run detection.
            %   All name-value options forwarded to BlinkAnalysis.compute():
            %     'method'      'mad' | 'mean'   (default: 'mad')
            %     'bandpass'    [lo hi] Hz        (default: [1 10])
            %     'multiplier'  scalar            (default: 5/6 per method)
            %     'minDist'     seconds           (default: 0.5/0.2 per method)
            obj.checkLoaded();
            obj.blink = BlinkAnalysis( ...
                obj.raw, obj.subject_results, obj.fs, ...
                obj.subjectID, obj.algorithmName);
            obj.blink.compute(varargin{:});
        end


        % ================================================================
        % computeSignalStats — delegates to SignalStatistics
        % ================================================================
        function computeSignalStats(obj)
            % COMPUTESIGNALSTATS - Instantiate SignalStatistics and compute.
            %   Segments: calibration, closed, open, all
            %   Window groups: trivial, nontrivial (from modifiedMask)
            %   Probe summary: trivialRate, meanNormR, meanRiemannDrift
            obj.checkLoaded();
            obj.signalStats = SignalStatistics( ...
                obj.raw, obj.subject_results, obj.fs, ...
                obj.subjectID, obj.algorithmName);
            obj.signalStats.compute();
        end

        % ================================================================
        % computeICA — delegates to ICAAnalysis
        % ================================================================
        function computeICA(obj, varargin)
            % COMPUTEICA - Instantiate ICAAnalysis and run IC-based evaluation.
            %
            %   Implements the Chang et al. (2020) framework: ICA decomposition
            %   of raw and cleaned signals, IC classification via ICLabel, and
            %   source power retention quantified per IC class.
            %
            %   Requires EEGLAB in the path (runica + iclabel).
            %   See ICAAnalysis for full option documentation.
            %
            %   Name-value options (all forwarded to ICAAnalysis.compute()):
            %     'segment'            'closed+open' (default) | 'closed' | 'open' | 'all'
            %     'maxSamples'         cap on samples per ICA run (default 60000)
            %     'preservationThresh' IC correlation threshold (default 0.8)
            %     'reuseRawW'          skip cleaned ICA, apply raw W only (default false)
            %     'chanlocs'           EEGLAB chanlocs struct for ICLabel (optional)
            %
            %   Example:
            %     ana.computeICA();
            %     ana.computeICA('segment','closed+open','maxSamples',30000);
            %     ana.ica.getSummary('powerRetained_brain')
            %     ana.ica.plotPreservation()
            %     ana.ica.plotPowerRetained()

            obj.checkLoaded();

            % Extract chanlocs from varargin if provided
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'chanlocs', [], @(x) isstruct(x)||isempty(x));
            parse(p, varargin{:});
            chanlocs = p.Results.chanlocs;

            % Remove 'chanlocs' from varargin so ICAAnalysis.compute() doesn't see it
            remaining = p.Unmatched;
            fwd = {};
            fn  = fieldnames(remaining);
            for i = 1:numel(fn)
                fwd{end+1} = fn{i};          %#ok<AGROW>
                fwd{end+1} = remaining.(fn{i}); %#ok<AGROW>
            end

            obj.ica = ICAAnalysis( ...
                obj.raw, obj.subject_results, obj.fs, ...
                obj.subjectID, obj.algorithmName, chanlocs);
            obj.ica.compute(fwd{:});
        end


        function vals = getSummary(obj, metricName, fieldName)
            % GETSUMMARY - Pull one scalar summary field across all combos.
            %   metricName : analyser name, e.g. 'blink'
            %   fieldName  : field in .summary, e.g. 'blinkReductionRatio_closed'
            switch metricName
                case 'blink'
                    if isempty(obj.blink)
                        error('ExperimentAnalysis:notComputed', ...
                            'Call computeBlinks() first.');
                    end
                    vals = obj.blink.getSummary(fieldName);
                case 'ica'
                    if isempty(obj.ica)
                        error('ExperimentAnalysis:notComputed', ...
                            'Call computeICA() first.');
                    end
                    vals = obj.ica.getSummary(fieldName);
                case 'signalStats'
                    if isempty(obj.signalStats)
                        error('ExperimentAnalysis:notComputed', ...
                            'Call computeSignalStats() first.');
                    end
                    % fieldName format: '<segment>.<signalType>.<field>'
                    % e.g. 'closed.clean.rms'  or  'probe.trivialRate'
                    parts = strsplit(fieldName, '.');
                    vals  = zeros(1, obj.nCombos);
                    for c = 1:obj.nCombos
                        r = obj.signalStats.results(c);
                        if numel(parts) == 2
                            vals(c) = r.(parts{1}).(parts{2});
                        elseif numel(parts) == 3
                            vals(c) = r.(parts{1}).(parts{2}).(parts{3});
                        end
                    end
                otherwise
                    error('ExperimentAnalysis:unknownMetric', ...
                        'Unknown metric "%s". Run the corresponding compute* method first.', ...
                        metricName);
            end
        end

        % ================================================================
        % paramTable — convenience passthrough
        % ================================================================
        function T = paramTable(obj)
            % PARAMTABLE - One row per combo with all scalar parameter values.
            obj.checkLoaded();
            if ~isempty(obj.blink)
                T = obj.blink.paramTable();
                return;
            end
            % Fallback: build directly from subject_results
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
                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

        % ================================================================
        % metrics — backward-compatible shim
        % ================================================================
        function m = metrics(obj, k)
            % METRICS - Backward-compatible access: ana.metrics(k).blink
            %   Preferred new access: ana.blink.results(k)
            if nargin < 2
                m = struct();
                for c = 1:obj.nCombos
                    if ~isempty(obj.blink) && numel(obj.blink.results) >= c
                        m(c).blink = obj.blink.results(c);
                    end
                end
            else
                m = struct();
                if ~isempty(obj.blink) && numel(obj.blink.results) >= k
                    m.blink = obj.blink.results(k);
                end
            end
        end

    end % public methods

    % ================================================================
    % Static factories
    % ================================================================
    methods (Static)

        function obj = fromSingleRun(resultsFolder, algorithmName, subjectID, fs, params)
            % FROMSINGLERUN - Single-run mode: reads one result .mat from disk.
            if ~isfolder(resultsFolder)
                error('ExperimentAnalysis:invalidFolder', ...
                    'resultsFolder does not exist: %s', resultsFolder);
            end
            if nargin < 5 || isempty(params), params = struct(); end

            obj                   = ExperimentAnalysis(resultsFolder, algorithmName, subjectID, fs);
            obj.accumulatedFolder = '';
            obj.resultsFolder     = resultsFolder;
            obj.mode              = 'single';
            obj.singleParams      = params;
        end

        function obj = fromOutputResults(data, asrObj, algorithmName, fs)
            % FROMOUTPUTRESULTS - Memory mode: no file I/O.
            obj               = ExperimentAnalysis(tempdir(), algorithmName, 0, fs);
            obj.mode          = 'memory';
            obj.algorithmName = algorithmName;

            obj.raw = struct( ...
                'calibration', data.calibrationSignal, ...
                'closed',      data.closedSignal,      ...
                'open',        data.openSignal);

            if isprop(asrObj, 'config') && ~isempty(asrObj.config)
                params = asrObj.config.algorithm.parameters;
            else
                params = struct();
            end

            sr = struct( ...
                'parameters',       params,                      ...
                'timestamp',        datetime('now'),             ...
                'cleanCalibration', data.cleanedCalibration,     ...
                'cleanClosed',      data.cleanedClosed,          ...
                'cleanOpen',        data.cleanedOpen,            ...
                'modifiedMask',     asrObj.modifiedMask,         ...
                'timeProbe',        [],                          ...
                'probeRaw',         [],                          ...
                'probeClean',       []);

            if isprop(asrObj, 'timeProbe') && ~isempty(asrObj.timeProbe)
                sr.timeProbe = asrObj.timeProbe;
            elseif isprop(asrObj, 'tProbe') && ~isempty(asrObj.tProbe)
                sr.timeProbe = asrObj.tProbe;
            end
            if isprop(asrObj, 'probeRaw'),   sr.probeRaw   = asrObj.probeRaw;   end
            if isprop(asrObj, 'probeClean'), sr.probeClean = asrObj.probeClean; end
            if isprop(asrObj, 'emaProbe'),   sr.probeRaw   = asrObj.emaProbe;   end

            obj.subject_results = sr;
            obj.nCombos         = 1;

            fprintf('Loaded (memory) %s -- ready.\n', algorithmName);
        end

    end % static methods

    % ================================================================
    % Private helpers
    % ================================================================
    methods (Access = private)

        function loadAccumulated(obj)
            fname = sprintf('S%d_%s_accumulated.mat', obj.subjectID, obj.algorithmName);
            fpath = fullfile(obj.accumulatedFolder, fname);
            if ~isfile(fpath)
                error('ExperimentAnalysis:notFound', ...
                    'Accumulated file not found: %s', fpath);
            end
            loaded              = load(fpath, 'raw', 'subject_results');
            obj.raw             = loaded.raw;
            obj.subject_results = loaded.subject_results;
            obj.nCombos         = numel(obj.subject_results);
            fprintf('Loaded (sweep) S%d %s -- %d combination(s).\n', ...
                obj.subjectID, obj.algorithmName, obj.nCombos);
        end

        function loadSingleRun(obj)
            paramStr = ExperimentAnalysis.buildParamString(obj.singleParams);
            fname    = sprintf('S%d_%s%s.mat', obj.subjectID, obj.algorithmName, paramStr);
            fpath    = fullfile(obj.resultsFolder, fname);
            if ~isfile(fpath)
                error('ExperimentAnalysis:notFound', ...
                    'Result file not found: %s', fpath);
            end
            loaded = load(fpath, 'results');
            r      = loaded.results;

            obj.raw = struct( ...
                'calibration', r.data.raw.cal,    ...
                'closed',      r.data.raw.closed, ...
                'open',        r.data.raw.open);

            sr = struct( ...
                'parameters',       r.config.algorithm.parameters, ...
                'timestamp',        r.timestamp,                   ...
                'cleanCalibration', r.data.cleaned.cal,            ...
                'cleanClosed',      r.data.cleaned.closed,         ...
                'cleanOpen',        r.data.cleaned.open,           ...
                'modifiedMask',     r.modifiedMask,                ...
                'timeProbe',        [],                            ...
                'probeRaw',         [],                            ...
                'probeClean',       []);

            if isfield(r, 'timeProbe'), sr.timeProbe = r.timeProbe; end
            if isfield(r, 'probes')
                if isfield(r.probes, 'raw'),   sr.probeRaw   = r.probes.raw;   end
                if isfield(r.probes, 'clean'), sr.probeClean = r.probes.clean; end
                if isfield(r.probes, 'ema'),   sr.probeRaw   = r.probes.ema;   end
            end

            obj.subject_results = sr;
            obj.nCombos         = 1;
            fprintf('Loaded (single) S%d %s -- params: %s\n', ...
                obj.subjectID, obj.algorithmName, paramStr);
        end

        function checkLoaded(obj)
            if isempty(obj.subject_results)
                error('ExperimentAnalysis:notLoaded', 'Call load() before computing metrics.');
            end
        end

    end % private methods

    methods (Static, Access = private)

        function str = buildParamString(params)
            fields = fieldnames(params);
            str    = "";
            for i = 1:numel(fields)
                val = params.(fields{i});
                if iscell(val), val = val{1}; end
                if isnumeric(val)
                    if isscalar(val)
                        valStr = num2str(val);
                    else
                        valStr = sprintf('%g-', val);
                        valStr(end) = [];
                    end
                else
                    valStr = char(string(val));
                end
                str = str + "_" + fields{i} + "_" + valStr;
            end
        end

    end % static private methods

end
