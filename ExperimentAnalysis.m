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
    %   MEMORY MODE — from pre-loaded data structs:
    %       ana = ExperimentAnalysis.fromData(raw, subject_results, 'vanilla-asr');
    %       ana = ExperimentAnalysis.fromData(raw, subject_results, 'vanilla-asr', 500, 3);
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
        spectral          % SpectralAnalysis instance
        quietRegion       % QuietRegionAnalysis instance
        temporal          % TemporalAnalysis instance
        ica               % ICAAnalysis instance

        % ---- Streaming state (populated by loadMeta) ----
        paramsList        % 1×nCombos cell array of param structs (no signal data)
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
            obj.blink.loaderFcn   = @(k) obj.loadCombo(k);
            obj.blink.unloaderFcn = @(k) obj.unloadCombo(k);
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
            obj.signalStats.loaderFcn   = @(k) obj.loadCombo(k);
            obj.signalStats.unloaderFcn = @(k) obj.unloadCombo(k);
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

        % ================================================================
        % computeSpectral — delegates to SpectralAnalysis
        % ================================================================
        function computeSpectral(obj, varargin)
            % COMPUTESPECTRAL - Instantiate SpectralAnalysis and run.
            %
            %   All name-value options forwarded to SpectralAnalysis.compute():
            %     'welchWinSec'    scalar > 0   Welch window length (sec)
            %                                   [] = MATLAB default  (default [])
            %     'welchOverlap'   0 ≤ x < 1    Overlap fraction
            %                                   [] = MATLAB default 50%  (default [])
            %     'freqResolution' scalar > 0   Freq-axis step Hz   (default 0.5)
            %
            %   Example:
            %     ana.computeSpectral();
            %     ana.computeSpectral('welchWinSec', 2, 'welchOverlap', 0.5);
            %     ana.spectral.getSummary('alpha.totalPowerAfter', 'closed')
            %     ana.spectral.plotPSD(1, 'closed')
            %     T = ana.spectral.paramTable();

            obj.checkLoaded();
            obj.spectral = SpectralAnalysis( ...
                obj.raw, obj.subject_results, obj.fs, ...
                obj.subjectID, obj.algorithmName);
            obj.spectral.loaderFcn   = @(k) obj.loadCombo(k);
            obj.spectral.unloaderFcn = @(k) obj.unloadCombo(k);
            obj.spectral.compute(varargin{:});
        end

        % ================================================================
        % computeQuietRegion — delegates to QuietRegionAnalysis
        % ================================================================
        function computeQuietRegion(obj, varargin)
            % COMPUTEQUIETREGION - Instantiate QuietRegionAnalysis and run.
            %   Requires computeBlinks() to have been called first.
            %
            %   All name-value options forwarded to QuietRegionAnalysis.compute():
            %     'quietWinSec'      scalar   Exclusion half-window around blinks (sec)
            %                                 default 0.25
            %     'maxDelaySec'      scalar   finddelay search range (sec)
            %                                 default 0.25
            %     'minQuietSamples'  scalar   Skip channel if fewer quiet samples
            %                                 default 50
            %
            %   Example:
            %     ana.computeBlinks();
            %     ana.computeQuietRegion();
            %     ana.computeQuietRegion('quietWinSec', 0.5);
            %     ana.quietRegion.getSummary('corr_median',   'closed')
            %     ana.quietRegion.getSummary('rre_median',    'closed')
            %     T = ana.quietRegion.paramTable();

            obj.checkLoaded();
            if isempty(obj.blink)
                error('ExperimentAnalysis:blinkRequired', ...
                    ['computeQuietRegion() requires blink data. ' ...
                     'Call computeBlinks() first.']);
            end
            obj.quietRegion = QuietRegionAnalysis( ...
                obj.raw, obj.subject_results, obj.fs, obj.blink, ...
                obj.subjectID, obj.algorithmName);
            obj.quietRegion.loaderFcn   = @(k) obj.loadCombo(k);
            obj.quietRegion.unloaderFcn = @(k) obj.unloadCombo(k);
            obj.quietRegion.compute(varargin{:});
        end

        % ================================================================
        % computeTemporal — delegates to TemporalAnalysis
        % ================================================================
        function computeTemporal(obj, varargin)
            % COMPUTETEMPORAL - Instantiate TemporalAnalysis and run.
            %
            %   All name-value options forwarded to TemporalAnalysis.compute():
            %     'windowSec'        scalar > 0  Window length (sec)      default 2
            %     'hopSec'           scalar > 0  Hop size (sec)           default 1
            %     'overlap'          0 ≤ x < 1   Overlap fraction         (overrides hopSec)
            %     'subspaceRank'     scalar       Top-k eigenvectors       default min(10,nCh)
            %     'channels'         vector       Channel subset           default all
            %     'compactionThresh' 0–1          Energy compaction target default 0.9
            %     'kRef'             scalar       Reference subspace rank  default subspaceRank
            %     'useRobustCov'     logical      Use block_geometric_median default true
            %     'blockSize'        scalar       Block size for robust cov  default 10
            %
            %   Example:
            %     ana.computeTemporal();
            %     ana.computeTemporal('windowSec', 2, 'overlap', 0.5, 'subspaceRank', 10);
            %     ana.temporal.getSummary('signal.rrmse_mean',        'closed')
            %     ana.temporal.getSummary('subspace.angle_mean_mean', 'closed')
            %     ana.temporal.plotTimeSeries(1, 'signal.rrmse', 'closed')
            %     T = ana.temporal.paramTable();

            obj.checkLoaded();
            obj.temporal = TemporalAnalysis( ...
                obj.raw, obj.subject_results, obj.fs, ...
                obj.subjectID, obj.algorithmName);
            obj.temporal.loaderFcn   = @(k) obj.loadCombo(k);
            obj.temporal.unloaderFcn = @(k) obj.unloadCombo(k);
            obj.temporal.compute(varargin{:});
        end


        % ================================================================
        % loadMeta — load raw signals + meta only (no combo data)
        % ================================================================
        function loadMeta(obj)
            % LOADMETA - Load raw signals and combo count/params without
            %   loading any cleaned signal data. Used for streaming workflows
            %   where combos are processed one at a time via loadCombo(k).
            %
            %   After loadMeta():
            %     obj.raw        -- populated
            %     obj.nCombos    -- populated
            %     obj.paramsList -- populated (cell array of param structs)
            %     obj.subject_results -- empty (not loaded yet)
            %
            %   Example:
            %     ana = ExperimentAnalysis(accumulatedDir, 'graph-asr', 6, 500);
            %     ana.loadMeta();
            %     for k = 1:ana.nCombos
            %         ana.loadCombo(k);
            %         ana.computeBlinks();
            %         % ...
            %     end

            rawFname  = sprintf('S%d_%s_raw.mat',  obj.subjectID, obj.algorithmName);
            metaFname = sprintf('S%d_%s_meta.mat', obj.subjectID, obj.algorithmName);
            rawFpath  = fullfile(obj.accumulatedFolder, rawFname);
            metaFpath = fullfile(obj.accumulatedFolder, metaFname);

            if ~isfile(rawFpath) || ~isfile(metaFpath)
                error('ExperimentAnalysis:noSplitLayout', ...
                    ['loadMeta() requires split layout files.\n' ...
                     'Expected:\n  %s\n  %s\n' ...
                     'Use load() for monolithic files.'], rawFpath, metaFpath);
            end

            loadedRaw      = load(rawFpath, 'raw');
            obj.raw        = loadedRaw.raw;

            loadedMeta     = load(metaFpath, 'nCombos', 'paramsList');
            obj.nCombos    = loadedMeta.nCombos;
            obj.paramsList = loadedMeta.paramsList;

            % Build a lightweight subject_results stub — params only, no signals.
            % This lets paramTable() and identity columns work before any combo loads.
            comboFields = {'parameters','timestamp', ...
                           'cleanCalibration','cleanClosed','cleanOpen', ...
                           'modifiedMask','timeProbe','probeRaw','probeClean'};
            prototype          = cell2struct(repmat({[]}, numel(comboFields), 1), comboFields);
            obj.subject_results = repmat(prototype, 1, obj.nCombos);
            for k = 1:obj.nCombos
                obj.subject_results(k).parameters = obj.paramsList{k};
            end

            fprintf('loadMeta: S%d %s -- %d combos ready for streaming.\n', ...
                obj.subjectID, obj.algorithmName, obj.nCombos);
        end

        % ================================================================
        % loadCombo — load one combo's cleaned signals into subject_results(k)
        % ================================================================
        function loadCombo(obj, k)
            % LOADCOMBO - Stream one combo file into subject_results(k).
            %   Replaces whatever was previously in subject_results(k).
            %   All other combos remain as they were (stub or previously loaded).
            %
            %   Requires loadMeta() to have been called first.
            %
            %   Inputs:
            %     k : combo index (1-based, 1 <= k <= obj.nCombos)
            %
            %   Example:
            %     ana.loadMeta();
            %     ana.loadCombo(3);
            %     % now only combo 3 has real signal data

            if isempty(obj.nCombos) || obj.nCombos == 0
                error('ExperimentAnalysis:notLoaded', ...
                    'Call loadMeta() before loadCombo().');
            end
            if k < 1 || k > obj.nCombos
                error('ExperimentAnalysis:outOfRange', ...
                    'k=%d out of range [1, %d].', k, obj.nCombos);
            end

            comboFname = sprintf('S%d_%s_combo_%d.mat', ...
                obj.subjectID, obj.algorithmName, k);
            comboFpath = fullfile(obj.accumulatedFolder, comboFname);

            if ~isfile(comboFpath)
                error('ExperimentAnalysis:comboNotFound', ...
                    'Combo file missing: %s', comboFpath);
            end

            loadedCombo = load(comboFpath, 'combo');
            c           = loadedCombo.combo;

            comboFields = {'parameters','timestamp', ...
                           'cleanCalibration','cleanClosed','cleanOpen', ...
                           'modifiedMask','timeProbe','probeRaw','probeClean'};
            for f = 1:numel(comboFields)
                fld = comboFields{f};
                if isfield(c, fld)
                    obj.subject_results(k).(fld) = c.(fld);
                end
            end
        end

        % ================================================================
        % unloadCombo — free signal data for combo k to recover memory
        % ================================================================
        function unloadCombo(obj, k)
            % UNLOADCOMBO - Clear cleaned signal data for combo k.
            %   Parameters stub is preserved so getSummary() identity columns work.
            %   Call after saving results for k to free memory before loading k+1.

            if k < 1 || k > obj.nCombos, return; end
            obj.subject_results(k).cleanCalibration = [];
            obj.subject_results(k).cleanClosed      = [];
            obj.subject_results(k).cleanOpen        = [];
            obj.subject_results(k).modifiedMask     = [];
            obj.subject_results(k).timeProbe        = [];
            obj.subject_results(k).probeRaw         = [];
            obj.subject_results(k).probeClean       = [];
        end


        function vals = getSummary(obj, metricName, fieldName, segment)
            % GETSUMMARY - Pull one scalar summary field across all combos.
            %
            %   metricName : analyser name
            %                'blink' | 'signalStats' | 'spectral' |
            %                'temporal' | 'quietRegion' | 'ica'
            %   fieldName  : field path — format depends on analyser (see each class)
            %   segment    : 'closed' | 'open' | 'all'  (required for spectral,
            %                temporal, quietRegion; ignored for blink, ica)
            %
            %   Examples:
            %     ana.getSummary('blink',       'blinkReductionRatio_closed')
            %     ana.getSummary('spectral',    'alpha.totalPowerAfter',    'closed')
            %     ana.getSummary('temporal',    'signal.rrmse_mean',        'closed')
            %     ana.getSummary('quietRegion', 'corr_median',              'closed')

            if nargin < 4, segment = 'closed'; end

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
                case 'spectral'
                    if isempty(obj.spectral)
                        error('ExperimentAnalysis:notComputed', ...
                            'Call computeSpectral() first.');
                    end
                    vals = obj.spectral.getSummary(fieldName, segment);
                case 'temporal'
                    if isempty(obj.temporal)
                        error('ExperimentAnalysis:notComputed', ...
                            'Call computeTemporal() first.');
                    end
                    vals = obj.temporal.getSummary(fieldName, segment);
                case 'quietRegion'
                    if isempty(obj.quietRegion)
                        error('ExperimentAnalysis:notComputed', ...
                            'Call computeQuietRegion() first.');
                    end
                    vals = obj.quietRegion.getSummary(fieldName, segment);
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

        function obj = fromData(raw, subject_results, algorithmName, fs, subjectID)
            % FROMDATA - Memory mode: construct directly from pre-loaded data.
            %
            %   obj = ExperimentAnalysis.fromData(raw, subject_results, algorithmName)
            %   obj = ExperimentAnalysis.fromData(raw, subject_results, algorithmName, fs)
            %   obj = ExperimentAnalysis.fromData(raw, subject_results, algorithmName, fs, subjectID)
            %
            %   Inputs:
            %     raw             - struct with fields .calibration, .closed, .open [C x T]
            %     subject_results - 1xN struct array of per-combo results
            %     algorithmName   - char, e.g. 'vanilla-asr'
            %     fs              - sampling frequency in Hz (default: 500)
            %     subjectID       - numeric subject identifier (default: 0)
            %
            %   This factory bypasses all file I/O. Useful when data is already
            %   in memory (e.g. from a batch view or external loader).
            if nargin < 5 || isempty(subjectID), subjectID = 0;   end
            if nargin < 4 || isempty(fs),        fs        = 500; end
            obj                   = ExperimentAnalysis(tempdir(), algorithmName, subjectID, fs);
            obj.mode              = 'memory';
            obj.raw               = raw;
            obj.subject_results   = subject_results;
            obj.nCombos           = numel(subject_results);
        end

    end % static methods

    % ================================================================
    % Private helpers
    % ================================================================
    methods (Access = private)

        function loadAccumulated(obj)
            % LOADACCUMULATED - Load sweep data from disk.
            %   Auto-detects split layout (_raw.mat + _meta.mat + _combo_k.mat)
            %   or falls back to the legacy monolithic _accumulated.mat.
            %
            %   Split layout (preferred):
            %     Reads raw and meta eagerly; assembles subject_results by
            %     streaming each _combo_k.mat in sequence. Memory footprint is
            %     identical to monolithic once all combos are loaded, but the
            %     files on disk stay small enough to avoid the 12 GB ceiling.
            %
            %   Monolithic layout (legacy):
            %     Loads raw + subject_results from a single _accumulated.mat.

            rawFname  = sprintf('S%d_%s_raw.mat',  obj.subjectID, obj.algorithmName);
            metaFname = sprintf('S%d_%s_meta.mat', obj.subjectID, obj.algorithmName);
            rawFpath  = fullfile(obj.accumulatedFolder, rawFname);
            metaFpath = fullfile(obj.accumulatedFolder, metaFname);

            if isfile(rawFpath) && isfile(metaFpath)
                % ── Split layout ──────────────────────────────────────
                obj.loadSplit(rawFpath, metaFpath);
            else
                % ── Monolithic layout (legacy) ─────────────────────────
                obj.loadMonolithic();
            end
        end

        function loadSplit(obj, rawFpath, metaFpath)
            % LOADSPLIT - Load split-layout files for this subject.
            %   Reads _raw.mat and _meta.mat eagerly, then streams each
            %   _combo_k.mat in sequence to assemble subject_results.
            %
            %   Combo files are expected at:
            %     <accumulatedFolder>/S<id>_<algo>_combo_<k>.mat
            %   Each contains a single variable 'combo' with fields:
            %     parameters, timestamp, cleanCalibration, cleanClosed,
            %     cleanOpen, modifiedMask, timeProbe, probeRaw, probeClean

            % ── Raw signals ───────────────────────────────────────────
            loadedRaw   = load(rawFpath, 'raw');
            obj.raw     = loadedRaw.raw;

            % ── Meta: nCombos + parameter list ────────────────────────
            loadedMeta  = load(metaFpath, 'nCombos', 'paramsList');
            nC          = loadedMeta.nCombos;

            % ── Combos: stream each file ──────────────────────────────
            comboFields = {'parameters','timestamp', ...
                           'cleanCalibration','cleanClosed','cleanOpen', ...
                           'modifiedMask','timeProbe','probeRaw','probeClean'};

            % Pre-allocate struct array using the field list
            prototype = cell2struct(repmat({[]}, numel(comboFields), 1), comboFields);
            sr        = repmat(prototype, 1, nC);

            for k = 1:nC
                comboFname = sprintf('S%d_%s_combo_%d.mat', ...
                    obj.subjectID, obj.algorithmName, k);
                comboFpath = fullfile(obj.accumulatedFolder, comboFname);

                if ~isfile(comboFpath)
                    error('ExperimentAnalysis:comboNotFound', ...
                        'Split combo file missing for S%d %s combo %d.\nExpected: %s', ...
                        obj.subjectID, obj.algorithmName, k, comboFpath);
                end

                loadedCombo = load(comboFpath, 'combo');
                c           = loadedCombo.combo;

                for f = 1:numel(comboFields)
                    fld = comboFields{f};
                    if isfield(c, fld)
                        sr(k).(fld) = c.(fld);
                    end
                end

                fprintf('  combo %d/%d loaded.\n', k, nC);
            end

            obj.subject_results = sr;
            obj.nCombos         = nC;
            fprintf('Loaded (sweep/split) S%d %s -- %d combination(s).\n', ...
                obj.subjectID, obj.algorithmName, nC);
        end

        function loadMonolithic(obj)
            % LOADMONOLITHIC - Load legacy single _accumulated.mat file.
            fname = sprintf('S%d_%s_accumulated.mat', obj.subjectID, obj.algorithmName);
            fpath = fullfile(obj.accumulatedFolder, fname);
            if ~isfile(fpath)
                error('ExperimentAnalysis:notFound', ...
                    ['No accumulated data found for S%d %s.\n' ...
                     'Checked split : %s_raw.mat / _meta.mat\n' ...
                     'Checked mono  : %s'], ...
                    obj.subjectID, obj.algorithmName, ...
                    sprintf('S%d_%s', obj.subjectID, obj.algorithmName), ...
                    fpath);
            end
            loaded              = load(fpath, 'raw', 'subject_results');
            obj.raw             = loaded.raw;
            obj.subject_results = loaded.subject_results;
            obj.nCombos         = numel(obj.subject_results);
            fprintf('Loaded (sweep/mono) S%d %s -- %d combination(s).\n', ...
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
