classdef Experimenting < handle
    % Experimenting - Orchestration engine for ASR benchmarking.
    %   Manages data flow between raw EEG datasets and the four ASR variants
    %   (originalASR, vanillaASR, emaASR, graphASR). Partitions data into
    %   functional segments (Calibration, Closed-Eyes, Open-Eyes) and captures
    %   telemetry via a shared timeProbe and per-variant diagnostic probes.

    properties (Access = public)
        config          % config_experiment object containing paths and params
        data            % Struct containing raw, filtered, and segmented signals
        asr             % Handle to polymorphic ASR object
        segmentIndices  % Tracking structure for temporal slice definitions
    end

    methods

        % ====================================================
        % Constructor & Lifecycle
        % ====================================================

        function obj = Experimenting(config)
            % EXPERIMENTING - Initialise the experiment runner.
            %   Input: config - A config_experiment instance.
            obj.config = config;

            obj.data = struct( ...
                'rawSignal',          [], ...
                'filteredSignal',     [], ...
                'calibrationSignal',  [], ...
                'closedSignal',       [], ...
                'openSignal',         [], ...
                'cleanedCalibration', [], ...
                'cleanedClosed',      [], ...
                'cleanedOpen',        []);

            obj.segmentIndices = struct('calibration', [], 'closed', [], 'open', []);
        end

        function obj = run(obj)
            % RUN - Execute the full experimental pipeline.
            fprintf('--- Starting Experiment: Subject %d [%s] ---\n', ...
                obj.config.subjectID, obj.config.algorithm.name);

            obj.loadData();
            obj.preprocess();
            obj.splitSegments();
            obj.calibrate();
            obj.process('all');
            obj.saveResults();

            fprintf('--- Experiment Complete ---\n');
        end

        % ====================================================
        % Data Ingestion & Preprocessing
        % ====================================================

        function obj = loadData(obj)
            % LOADDATA - Retrieve signal and metadata from the dataset directory.
            [signal, chanlabel, latency, chanlocs] = ...
                loadSubjectData(obj.config.subjectID, obj.config.paths.dataset);

            obj.data.rawSignal        = double(signal);
            obj.data.channelName      = chanlabel;
            obj.data.channelLocations = chanlocs;
            obj.data.timeEvent        = latency;
        end

        function obj = preprocess(obj)
            % PREPROCESS - Apply zero-phase filtering and average re-referencing.
            fs = obj.config.signal.fs;
            bp = obj.config.signal.bandpass;
            ln = obj.config.signal.lineNoise;

            x = obj.data.rawSignal;
            x = x - mean(x, 2);        % Common Average Reference (CAR)
            x = obj.applyBandpass(x, fs, bp);
            x = obj.applyNotch(x, fs, ln);

            obj.data.filteredSignal = x;
        end

        % ====================================================
        % Segmentation & Calibration
        % ====================================================

        function obj = splitSegments(obj)
            % SPLITSEGMENTS - Partition data into calibration and test sets.
            x         = obj.data.filteredSignal;
            fs        = obj.config.signal.fs;
            params    = obj.config.algorithm.parameters;
            timeEvent = obj.data.timeEvent;

            [cal, closed, open, padded, idx] = obj.prepareASRData( ...
                x, fs, params.calibWindSec, params.lookahead, timeEvent);

            obj.data.calibrationSignal  = cal;
            obj.data.closedSignal       = closed;
            obj.data.openSignal         = open;
            obj.data.contaminatedSignal = padded;
            obj.data.indices            = idx;
        end

        function obj = calibrate(obj)
            % CALIBRATE - Instantiate and train the specific ASR variant.
            fs      = obj.config.signal.fs;
            params  = obj.config.algorithm.parameters;
            calEEG  = obj.data.calibrationSignal;
            algName = lower(obj.config.algorithm.name);

            % FIX #1: was asrTimeProbe() — class is timeProbe
            % FIX #4: single shared probe instance passed to ALL variants
            tp = timeProbe();

            fprintf('--- [Calibration] Initializing: %s ---\n', algName);

            switch algName
                case {'original-asr', 'o-asr'}
                    % FIX #4: originalASR constructor updated to accept timeProbe
                    obj.asr = originalASR(calEEG, fs, params.cutoff, tp);

                case {'vanilla-asr', 'v-asr'}
                    obj.asr = vanillaASR(calEEG, fs, params.cutoff, params.blocksize, tp);

                case {'hmo_asr', 'hmo-asr'}
                    obj.asr = hmoASR_working(fs, params);

                case {'ema_asr', 'ema-asr'}
                    % FIX #4: pass shared tp as third argument
                    obj.asr = emaASR(fs, params, tp);

                case {'e_asr', 'embeddedasr', 'e-asr'}
                    obj.asr = embeddedASR(fs, params, tp);

                case {'r_asr', 'riemann-asr', 'r-asr'}
                    obj.asr = riemannASR(fs, params, tp);

                case {'g_asr', 'graph-asr', 'g-asr-f', 'g-asr', 'g_asr-frank'}
                    % FIX #4: pass shared tp as third argument
                    obj.asr = graphASR(fs, params, tp);

                case {'g_asr-b', 'graph-asr-base', 'g-asr-b', 'g-asr-base', 'g_asr-base'}
                    obj.asr = graphASR_base(fs, params);

                otherwise
                    error('Experimenting:UnknownAlgorithm', ...
                        'Algorithm "%s" is not supported.', algName);
            end

            obj.asr.calibrate(calEEG);
        end

        % ====================================================
        % Processing Loop
        % ====================================================

        function obj = process(obj, mode)
            % PROCESS - Apply subspace reconstruction to selected segments.
            if nargin < 2, mode = 'all'; end

            % FIX #10: guard against calling before calibrate()
            if isempty(obj.asr)
                error('Experimenting:notCalibrated', ...
                    'Call calibrate() before process().');
            end

            switch lower(mode)
                case 'calibration'
                    fprintf('--- [Processing] Cleaning: Calibration Signal ---\n');
                    obj.data.cleanedCalibration = obj.asr.process(obj.data.calibrationSignal);

                case 'closed'
                    fprintf('--- [Processing] Cleaning: Eyes Closed Signal ---\n');
                    obj.data.cleanedClosed = obj.asr.process(obj.data.closedSignal);

                case 'open'
                    fprintf('--- [Processing] Cleaning: Eyes Open Signal ---\n');
                    obj.data.cleanedOpen = obj.asr.process(obj.data.openSignal);

                case 'all'
                    % Sequential order is critical for adaptive/EMA state tracking
                    fprintf('--- [Processing] Cleaning: Calibration Signal ---\n');
                    obj.data.cleanedCalibration = obj.asr.process(obj.data.calibrationSignal);
                    fprintf('--- [Processing] Cleaning: Eyes Closed Signal ---\n');
                    obj.data.cleanedClosed      = obj.asr.process(obj.data.closedSignal);
                    fprintf('--- [Processing] Cleaning: Eyes Open Signal ---\n');
                    obj.data.cleanedOpen        = obj.asr.process(obj.data.openSignal);
                    fprintf('--- [Processing] Done ---\n');

                otherwise
                    error('Experimenting:unknownMode', ...
                        'Unknown process mode "%s". Use: calibration | closed | open | all.', mode);
            end
        end

        function [data, asrObj] = outputResults(obj)
            data   = obj.data;
            asrObj = obj.asr;
        end

        % ====================================================
        % Results & IO
        % ====================================================

        function savePath = saveResults(obj)
            % SAVERESULTS - Export experiment state and cleaned data to .mat.

            % FIX #11: guard against saving before processing is complete
            if isempty(obj.data.cleanedCalibration) || ...
               isempty(obj.data.cleanedClosed)      || ...
               isempty(obj.data.cleanedOpen)
                error('Experimenting:saveResults', ...
                    'Cleaned data is empty. Run process() before saveResults().');
            end

            subjID = obj.config.subjectID;
            alg    = obj.config.algorithm.name;
            params = obj.config.algorithm.parameters;

            paramStr = obj.buildParamString(params);
            fname    = sprintf('S%d_%s%s.mat', subjID, alg, paramStr);
            savePath = fullfile(obj.config.paths.results, fname);

            results.config       = obj.config;
            results.data.cleaned = struct( ...
                'cal',    obj.data.cleanedCalibration, ...
                'closed', obj.data.cleanedClosed,      ...
                'open',   obj.data.cleanedOpen);

            if obj.config.sweepMode
                % SWEEP MODE: raw signals are identical across all combos
                % for a given subject — store once separately, skip here.
                % accumulateResults reads the raw file when accumulating.
                obj.saveRawOnce();
            else
                % SINGLE RUN MODE: save raw inline so the file is self-contained.
                results.data.raw = struct( ...
                    'cal',    obj.data.calibrationSignal, ...
                    'closed', obj.data.closedSignal,      ...
                    'open',   obj.data.openSignal);
            end
            results.modifiedMask = obj.asr.modifiedMask;
            results.timestamp    = datetime('now');

            % FIX #2/#5: read timing probe via helper that handles both
            % property names (tProbe for emaASR/graphASR, timeProbe for
            % vanillaASR/originalASR) rather than hard-coding one name.
            results.timeProbe = obj.getProbeFromASR();

            % FIX #6: save variant-specific diagnostic probes
            % vanillaASR/originalASR → probeRaw, probeClean
            % emaASR                 → emaProbe
            % graphASR               → probeRaw, probeClean
            if isprop(obj.asr, 'probeRaw')
                results.probes = struct('raw', obj.asr.probeRaw, 'clean', obj.asr.probeClean);
            elseif isprop(obj.asr, 'emaProbe')
                results.probes = struct('ema', obj.asr.emaProbe);
            end

            save(savePath, 'results', '-v7.3');
            fprintf('--- [Results] Saved to: %s ---\n', savePath);
        end

    end % end public methods

    methods (Access = private)

        % ====================================================
        % Signal Processing Utilities
        % ====================================================

        function [cal, closed, open, padded, idx] = prepareASRData(~, x, fs, calWindow, lookahead, timeEvent)
            % Pad signal with zeros to accommodate lookahead buffer
            lookaheadSamples = round(lookahead * fs);
            padded = [x, zeros(size(x,1), lookaheadSamples)];
            L      = size(padded, 2);

            % Robust handling of calWindow as numeric array or cell
            if iscell(calWindow), cw = calWindow{1}; else, cw = calWindow; end

            cStart = max(1, round(cw(1) * fs));
            cEnd   = min(L, round(cw(2) * fs));

            % FIX #8: guard against calibration window overlapping the event marker
            if cEnd >= timeEvent
                error('Experimenting:segmentOverlap', ...
                    ['Calibration window end (%d samples) meets or exceeds timeEvent (%d). ' ...
                     'Reduce calibWindSec or check timeEvent.'], cEnd, timeEvent);
            end

            cal    = padded(:, cStart:cEnd);
            closed = padded(:, (cEnd+1):(timeEvent-1));
            open   = padded(:, timeEvent:L);

            idx = struct( ...
                'calibration', [cStart,    cEnd        ], ...
                'closed',      [cEnd+1,    timeEvent-1 ], ...
                'open',        [timeEvent, L           ]);
        end

        function y = applyBandpass(~, x, fs, bp)
            nyq = fs / 2;
            b   = fir1(128, bp / nyq, 'bandpass');
            y   = filtfilt(b, 1, x')';
        end

        function y = applyNotch(~, x, fs, notch)
            nyq = fs / 2;
            b   = fir1(128, [(notch-1)/nyq, (notch+1)/nyq], 'stop');
            y   = filtfilt(b, 1, x')';
        end

        function saveRawOnce(obj)
            % SAVERAWONCE - Write raw signals for this subject once.
            %   File: S<id>_<alg>_raw.mat in the results folder.
            %   Skipped silently if already written by a previous combo run.
            subjID = obj.config.subjectID;
            alg    = obj.config.algorithm.name;
            fname  = sprintf('S%d_%s_raw.mat', subjID, alg);
            fpath  = fullfile(obj.config.paths.results, fname);

            if isfile(fpath), return; end   % already written — skip

            raw = struct( ...
                'cal',    obj.data.calibrationSignal, ...
                'closed', obj.data.closedSignal,      ...
                'open',   obj.data.openSignal);       %#ok<NASGU>

            save(fpath, 'raw', '-v7.3');
        end

        function probe = getProbeFromASR(obj)
            % FIX #2/#5: unified accessor for the timing probe regardless of
            % which property name the ASR variant uses.
            if isprop(obj.asr, 'tProbe') && ~isempty(obj.asr.tProbe)
                probe = obj.asr.tProbe;
            elseif isprop(obj.asr, 'timeProbe') && ~isempty(obj.asr.timeProbe)
                probe = obj.asr.timeProbe;
            else
                probe = [];
                warning('Experimenting:noProbe', ...
                    'ASR variant has no timing probe property (tProbe or timeProbe).');
            end
        end

        function str = buildParamString(~, params)
            % FIX #9: normalised indentation to match rest of class
            fields = fieldnames(params);
            str    = "";

            for i = 1:length(fields)
                val = params.(fields{i});

                if iscell(val)
                    val = val{1};
                end

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

    end % end private methods

end
