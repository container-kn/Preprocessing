classdef accumulateResults < handle
    % accumulateResults - Accumulate ASR sweep results per subject.
    %   Scans a results folder for all .mat files matching a subject list
    %   and algorithm name. For each subject, collects all parameter
    %   combinations and saves one .mat per subject containing:
    %       raw             - struct with .calibration/.closed/.open (stored ONCE)
    %       subject_results - 1xN struct array, one entry per combo, containing
    %                         cleaned signals + diagnostics only
    %
    %   Usage:
    %       acc = accumulateResults(resultsFolder, 'vanilla-asr', [1 2 3 4 5]);
    %       acc.run();
    %
    %       R   = acc.load(3);      % subject_results for subject 3
    %       raw = acc.loadRaw(3);   % raw signals for subject 3
    %
    %       acc.cleanup();          % delete individual run files (all subjects)
    %       acc.cleanup(3);         % delete individual run files (subject 3 only)

    properties (Access = public)
        baseFolder      % Path to folder containing individual run .mat files
        algorithmName   % Algorithm name string to match in filenames
        subjects        % Numeric row vector of subject IDs to process
        outputFolder    % Where per-subject accumulated files are written
        files           % Cell array of all matched file paths
        subjectMap      % containers.Map: subjectID -> {file paths}
        reportEvery     % Print progress every N files parsed (default 5)
    end

    methods

        % ================================================================
        % Constructor
        % ================================================================
        function obj = accumulateResults(baseFolder, algorithmName, subjects)
            if ~isfolder(baseFolder)
                error('accumulateResults:invalidFolder', ...
                    'baseFolder does not exist: %s', baseFolder);
            end
            obj.baseFolder    = baseFolder;
            obj.algorithmName = algorithmName;
            obj.subjects      = subjects(:)';
            obj.outputFolder  = fullfile(baseFolder, 'accumulated');
            obj.files         = {};
            obj.subjectMap    = containers.Map('KeyType','double','ValueType','any');
            obj.reportEvery   = 5;
        end

        % ================================================================
        % Locate and index all matching files grouped by subject
        % ================================================================
        function locateFiles(obj)
            d = dir(fullfile(obj.baseFolder, '*.mat'));

            for i = 1:numel(d)
                fname  = d(i).name;

                % FIX #2: skip raw-signal sidecar files written by saveRawOnce()
                if contains(fname, '_raw.mat'), continue; end

                tokens = regexp(fname, 'S(\d+)_([^_]+)', 'tokens');
                if isempty(tokens), continue; end

                subjectID = str2double(tokens{1}{1});
                algName   = tokens{1}{2};

                if ~ismember(subjectID, obj.subjects), continue; end
                if ~strcmpi(algName, obj.algorithmName), continue; end

                fpath            = fullfile(d(i).folder, fname);
                obj.files{end+1} = fpath;

                if obj.subjectMap.isKey(subjectID)
                    existing        = obj.subjectMap(subjectID);
                    existing{end+1} = fpath;
                    obj.subjectMap(subjectID) = existing;
                else
                    obj.subjectMap(subjectID) = {fpath};
                end
            end

            if isempty(obj.files)
                warning('accumulateResults:noFiles', ...
                    'No files found in %s for algorithm "%s", subjects [%s].', ...
                    obj.baseFolder, obj.algorithmName, num2str(obj.subjects));
                return;
            end

            fprintf('Located %d files across %d subjects.\n', ...
                numel(obj.files), obj.subjectMap.Count);
        end

        % ================================================================
        % Parse one run .mat file
        % ================================================================
        function entry = parseFile(~, filePath)
            S = load(filePath);
            if ~isfield(S, 'results')
                error('accumulateResults:missingField', ...
                    'File %s has no results struct.', filePath);
            end
            r = S.results;

            % FIX #1: raw fields only present in single-run mode files.
            % In sweep mode, raw signals live in the separate _raw.mat sidecar
            % written by saveRawOnce() -- read as empty here, filled later.
            if isfield(r.data, 'raw')
                rawCal    = r.data.raw.cal;
                rawClosed = r.data.raw.closed;
                rawOpen   = r.data.raw.open;
            else
                rawCal    = [];
                rawClosed = [];
                rawOpen   = [];
            end

            entry = struct( ...
                'subjectID',        r.config.subjectID,             ...
                'algorithm',        r.config.algorithm.name,        ...
                'parameters',       r.config.algorithm.parameters,  ...
                'timestamp',        r.timestamp,                    ...
                'rawCalibration',   rawCal,                         ...
                'rawClosed',        rawClosed,                      ...
                'rawOpen',          rawOpen,                        ...
                'cleanCalibration', r.data.cleaned.cal,             ...
                'cleanClosed',      r.data.cleaned.closed,          ...
                'cleanOpen',        r.data.cleaned.open,            ...
                'modifiedMask',     r.modifiedMask,                 ...
                'timeProbe',        [],                             ...
                'probeRaw',         [],                             ...
                'probeClean',       [],                             ...
                'config',           r.config);

            if isfield(r, 'timeProbe')
                entry.timeProbe = r.timeProbe;
            end
            if isfield(r, 'probes')
                if isfield(r.probes, 'raw'),   entry.probeRaw   = r.probes.raw;   end
                if isfield(r.probes, 'clean'), entry.probeClean = r.probes.clean; end
                if isfield(r.probes, 'ema'),   entry.probeRaw   = r.probes.ema;   end
            end
        end

        % ================================================================
        % Accumulate all combos for one subject and save
        % ================================================================
        function accumulateSubject(obj, subjectID, fileOffset, totalFiles, tStart)
            % fileOffset  : number of files already parsed before this subject
            % totalFiles  : total files across all subjects (for ETA)
            % tStart      : tic value from run() for elapsed time
            if nargin < 3, fileOffset  = 0;   end
            if nargin < 4, totalFiles  = [];   end
            if nargin < 5, tStart      = [];   end

            if ~obj.subjectMap.isKey(subjectID)
                warning('accumulateResults:unknownSubject', ...
                    'No files found for subject %d.', subjectID);
                return;
            end

            filePaths = obj.subjectMap(subjectID);
            N         = numel(filePaths);
            entries   = cell(N, 1);

            fprintf('  Subject %d -- %d combinations\n', subjectID, N);

            for i = 1:N
                try
                    entries{i} = obj.parseFile(filePaths{i});
                catch e
                    warning('accumulateResults:parseFailed', ...
                        'Skipping %s: %s', filePaths{i}, e.message);
                end

                % Progress + ETA every reportEvery files
                filesProcessedGlobal = fileOffset + i;
                if ~isempty(totalFiles) && ~isempty(tStart)
                    if mod(filesProcessedGlobal, obj.reportEvery) == 0 || ...
                            filesProcessedGlobal == totalFiles
                        obj.printProgress( ...
                            filesProcessedGlobal, totalFiles, ...
                            subjectID, i, N, tStart);
                    end
                end
            end

            valid   = ~cellfun(@isempty, entries);
            entries = entries(valid);
            nValid  = numel(entries);

            fprintf('    %d/%d parsed OK.\n', nValid, N);

            % Raw signals: read from the dedicated per-subject raw file
            % written once by Experimenting.saveRawOnce() during the sweep.
            rawFname = sprintf('S%d_%s_raw.mat', subjectID, obj.algorithmName);
            rawFpath = fullfile(obj.baseFolder, rawFname);
            if isfile(rawFpath)
                loaded = load(rawFpath, 'raw');
                raw    = struct( ...
                    'calibration', loaded.raw.cal,    ...
                    'closed',      loaded.raw.closed, ...
                    'open',        loaded.raw.open);
            else
                % Fallback: extract from first entry (sweep run with old Experimenting)
                warning('accumulateResults:noRawFile', ...
                    'Raw file not found (%s). Extracting from first combo entry.', rawFname);
                raw = struct( ...
                    'calibration', entries{1}.rawCalibration, ...
                    'closed',      entries{1}.rawClosed,      ...
                    'open',        entries{1}.rawOpen);
            end

            % Per combo: cleaned signals + diagnostics only
            subject_results = struct( ...
                'parameters',       [], ...
                'timestamp',        [], ...
                'cleanCalibration', [], ...
                'cleanClosed',      [], ...
                'cleanOpen',        [], ...
                'modifiedMask',     [], ...
                'timeProbe',        [], ...
                'probeRaw',         [], ...
                'probeClean',       []);

            subject_results = repmat(subject_results, 1, nValid);

            comboFields = {'parameters','timestamp', ...
                           'cleanCalibration','cleanClosed','cleanOpen', ...
                           'modifiedMask','timeProbe','probeRaw','probeClean'};

            for i = 1:nValid
                for f = 1:numel(comboFields)
                    subject_results(i).(comboFields{f}) = ...
                        entries{i}.(comboFields{f});
                end
            end

            if ~exist(obj.outputFolder, 'dir'), mkdir(obj.outputFolder); end

            fname = sprintf('S%d_%s_accumulated.mat', subjectID, obj.algorithmName);
            fpath = fullfile(obj.outputFolder, fname);
            save(fpath, 'raw', 'subject_results', '-v7.3');

            fprintf('    Saved: %s\n', fname);
        end

        % ================================================================
        % Run -- accumulate and save for all subjects
        % ================================================================
        function run(obj)
            obj.locateFiles();
            if isempty(obj.files), return; end

            if ~exist(obj.outputFolder, 'dir'), mkdir(obj.outputFolder); end

            % Total file count across all subjects for ETA
            totalFiles = numel(obj.files);

            fprintf('\n=== Accumulating: %s ===\n', obj.algorithmName);
            fprintf('  Total files   : %d\n', totalFiles);
            fprintf('  Subjects      : %d\n', numel(obj.subjects));
            fprintf('  Output folder : %s\n', obj.outputFolder);
            fprintf('  Progress every: %d files\n\n', obj.reportEvery);

            tStart = tic;

            fileOffset = 0;
            for s = 1:numel(obj.subjects)
                sid = obj.subjects(s);
                obj.accumulateSubject(sid, fileOffset, totalFiles, tStart);
                if obj.subjectMap.isKey(sid)
                    fileOffset = fileOffset + numel(obj.subjectMap(sid));
                end
            end

            elapsed = toc(tStart);
            fprintf('\nDone. %d subject files written in %.1fs.\n', ...
                numel(obj.subjects), elapsed);
        end

        % ================================================================
        % Load accumulated subject_results from disk
        % ================================================================
        function subject_results = load(obj, subjectID)
            % LOAD - Returns 1xN struct array, one entry per parameter combo.
            %   subject_results(k).parameters       -- param combo k
            %   subject_results(k).cleanOpen         -- [C x T] cleaned signal
            %   subject_results(k).modifiedMask      -- reconstruction mask
            %   subject_results(k).timeProbe         -- timing data
            fname = sprintf('S%d_%s_accumulated.mat', subjectID, obj.algorithmName);
            fpath = fullfile(obj.outputFolder, fname);
            if ~isfile(fpath)
                error('accumulateResults:notFound', ...
                    'No accumulated file for subject %d.\nExpected: %s\nRun run() first.', ...
                    subjectID, fpath);
            end
            loaded          = load(fpath, 'subject_results');
            subject_results = loaded.subject_results;
        end

        % ================================================================
        % Load raw signals for one subject
        % ================================================================
        function raw = loadRaw(obj, subjectID)
            % LOADRAW - Returns raw struct with .calibration/.closed/.open
            fname = sprintf('S%d_%s_accumulated.mat', subjectID, obj.algorithmName);
            fpath = fullfile(obj.outputFolder, fname);
            if ~isfile(fpath)
                error('accumulateResults:notFound', ...
                    'No accumulated file for subject %d. Run run() first.', subjectID);
            end
            loaded = load(fpath, 'raw');
            raw    = loaded.raw;
        end

        % ================================================================
        % Delete individual run files after accumulation
        % ================================================================
        function cleanup(obj, subjectID)
            % CLEANUP - Delete individual per-run .mat files from baseFolder.
            %   Checks that the accumulated file exists before deleting anything.
            %
            %   cleanup()          -- delete files for ALL subjects
            %   cleanup(subjectID) -- delete files for one subject only

            if nargin < 2
                targetIDs = obj.subjects;
            else
                targetIDs = subjectID;
            end

            nDeleted = 0;

            for s = 1:numel(targetIDs)
                sid = targetIDs(s);

                % Safety: confirm accumulated file exists before deleting
                fname = sprintf('S%d_%s_accumulated.mat', sid, obj.algorithmName);
                fpath = fullfile(obj.outputFolder, fname);
                if ~isfile(fpath)
                    warning('accumulateResults:cleanup', ...
                        'Accumulated file for subject %d not found -- skipping.\n  Expected: %s', ...
                        sid, fpath);
                    continue;
                end

                if ~obj.subjectMap.isKey(sid)
                    fprintf('  Subject %d -- no files indexed, skipping.\n', sid);
                    continue;
                end

                filePaths = obj.subjectMap(sid);
                fprintf('  Subject %d -- deleting %d files...\n', sid, numel(filePaths));

                for i = 1:numel(filePaths)
                    if isfile(filePaths{i})
                        delete(filePaths{i});
                        nDeleted = nDeleted + 1;
                    end
                end
            end

            fprintf('Cleanup done. %d files deleted.\n', nDeleted);
        end

    end % public methods

    methods (Access = private)

        function printProgress(obj, filesGlobal, totalFiles, subjectID, comboIdx, totalCombos, tStart)
            % PRINTPROGRESS - Print [X/N] status line with ETA.
            elapsed  = toc(tStart);
            rate     = filesGlobal / elapsed;           % files per second
            remaining = (totalFiles - filesGlobal) / max(rate, 1e-6);

            fprintf('  [%d/%d] S%d combo %d/%d | elapsed %s | ETA %s\n', ...
                filesGlobal, totalFiles, ...
                subjectID, comboIdx, totalCombos, ...
                obj.formatDuration(elapsed), ...
                obj.formatDuration(remaining));
        end

    end % private methods

    methods (Static, Access = private)

        function str = formatDuration(seconds)
            % FORMATDURATION - Convert seconds to mm:ss string.
            m   = floor(seconds / 60);
            s   = floor(mod(seconds, 60));
            str = sprintf('%02d:%02d', m, s);
        end

    end

end
