classdef accumulateResults < handle
    % accumulateResults - Accumulate ASR sweep results per subject.
    %   Scans a results folder for all .mat files matching a subject list
    %   and algorithm name. Supports two output layouts:
    %
    %   MONOLITHIC layout (splitFiles = false, legacy default):
    %       S<id>_<algo>_accumulated.mat  -- raw + all subject_results in one file
    %
    %   SPLIT layout (splitFiles = true, recommended for large sweeps):
    %       S<id>_<algo>_raw.mat          -- raw signals only
    %       S<id>_<algo>_meta.mat         -- nCombos + parameters for each combo
    %       S<id>_<algo>_combo_<k>.mat    -- cleaned signals + diagnostics for combo k
    %
    %   The split layout avoids 12+ GB monolithic files. ExperimentAnalysis
    %   detects the format automatically and streams combos on demand.
    %
    %   Usage:
    %       acc = accumulateResults(resultsFolder, 'vanilla-asr', [1 2 3 4 5]);
    %       acc.splitFiles = true;   % enable split layout (recommended)
    %       acc.run();               % first-time full accumulation
    %
    %       % --- Incremental workflow (run new combos, then patch) ---
    %       acc.patch();             % append only new combo files (skips duplicates)
    %       acc.patchCleanup();      % delete individual run files that were patched in
    %       acc.listCombos(3);       % print all accumulated combos for subject 3
    %
    %       R   = acc.load(3);           % subject_results for subject 3
    %       raw = acc.loadRaw(3);        % raw signals for subject 3
    %       sr  = acc.ensureComboLoaded(3, k);  % stream one combo (split layout)
    %
    %       acc.resplit(3);       % migrate monolithic file to split layout
    %       acc.resplit();        % migrate all subjects
    %
    %       acc.cleanup();        % delete individual run files (all subjects)
    %       acc.cleanup(3);       % delete individual run files (subject 3 only)

    properties (Access = public)
        baseFolder      % Path to folder containing individual run .mat files
        algorithmName   % Algorithm name string to match in filenames
        subjects        % Numeric row vector of subject IDs to process
        outputFolder    % Where per-subject accumulated files are written
        files           % Cell array of all matched file paths
        subjectMap      % containers.Map: subjectID -> {file paths}
        reportEvery     % Print progress every N files parsed (default 5)
        splitFiles      % true = write split layout; false = monolithic (default false)
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
            obj.splitFiles    = false;
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

            if obj.splitFiles
                obj.saveSplit(subjectID, raw, subject_results, nValid);
            else
                fname = sprintf('S%d_%s_accumulated.mat', subjectID, obj.algorithmName);
                fpath = fullfile(obj.outputFolder, fname);
                save(fpath, 'raw', 'subject_results', '-v7.3');
                fprintf('    Saved: %s\n', fname);
            end
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
        % ensureComboLoaded — stream one combo from split layout
        % ================================================================
        function sr = ensureComboLoaded(obj, subjectID, k)
            % ENSURECOMBOLOADED - Load combo k for subjectID from split layout.
            %   Returns the subject_results struct for combo k only.
            %   Raises an error if the split layout files are not present.
            %
            %   Inputs:
            %     subjectID : scalar integer subject ID
            %     k         : combo index (1-based)
            %
            %   Output:
            %     sr : struct with fields: parameters, timestamp,
            %          cleanCalibration, cleanClosed, cleanOpen,
            %          modifiedMask, timeProbe, probeRaw, probeClean
            %
            %   Example:
            %     sr3 = acc.ensureComboLoaded(5, 3);
            %     sr3.cleanClosed   % [C x T] cleaned closed-eyes signal

            fname = sprintf('S%d_%s_combo_%d.mat', subjectID, obj.algorithmName, k);
            fpath = fullfile(obj.outputFolder, fname);
            if ~isfile(fpath)
                error('accumulateResults:comboNotFound', ...
                    'Split combo file not found for S%d combo %d.\nExpected: %s\nRun run() with splitFiles=true first.', ...
                    subjectID, k, fpath);
            end
            loaded = load(fpath, 'combo');
            sr     = loaded.combo;
        end

        % ================================================================
        % patch — incrementally append new combo files to split layout
        % ================================================================
        function nAdded = patch(obj, subjectID)
            % PATCH - Scan for new individual combo files and append them to
            %   the existing split layout without touching combos already stored.
            %
            %   For each subject, this method:
            %     1. Reads the existing _meta.mat to get already-known params.
            %     2. Scans baseFolder for individual combo .mat files.
            %     3. Skips any file whose parameter set already exists in meta.
            %     4. Appends each new combo as _combo_<k+1>.mat, <k+2>.mat...
            %     5. Rewrites _meta.mat with the updated nCombos + paramsList.
            %
            %   If no split layout exists yet for a subject, falls back to a
            %   full run() for that subject first.
            %
            %   patch()            -- patch all subjects in obj.subjects
            %   patch(subjectID)   -- patch one subject only
            %
            %   Returns nAdded : total number of new combos written.
            %
            %   Example:
            %     acc = accumulateResults(resultsFolder, 'graph-asr', [6]);
            %     acc.splitFiles = true;
            %     acc.patch();     % run after any new sweep batch

            if nargin < 2
                targetIDs = obj.subjects;
            else
                targetIDs = subjectID(:)';
            end

            if ~exist(obj.outputFolder, 'dir'), mkdir(obj.outputFolder); end

            nAdded = 0;
            obj.locateFiles();

            for s = 1:numel(targetIDs)
                sid = targetIDs(s);
                fprintf('\n--- patch S%d %s ---\n', sid, obj.algorithmName);

                metaFname = sprintf('S%d_%s_meta.mat', sid, obj.algorithmName);
                metaFpath = fullfile(obj.outputFolder, metaFname);

                % ── Load existing meta (or bootstrap empty state) ──────
                wasBootstrapped = false;
                if isfile(metaFpath)
                    loaded         = load(metaFpath, 'nCombos', 'paramsList');
                    existingN      = loaded.nCombos;
                    existingParams = loaded.paramsList;   % 1×K cell of param structs
                    fprintf('  Existing combos : %d\n', existingN);
                else
                    % No split layout yet — full first-time build from all indexed files
                    fprintf('  No existing split layout found -- running full accumulation first.\n');
                    obj.splitFiles = true;
                    if obj.subjectMap.isKey(sid)
                        obj.accumulateSubject(sid);
                        loaded         = load(metaFpath, 'nCombos', 'paramsList');
                        existingN      = loaded.nCombos;
                        existingParams = loaded.paramsList;
                        fprintf('  Bootstrap done: %d combos written.\n', existingN);
                        wasBootstrapped = true;
                    else
                        warning('accumulateResults:patch', ...
                            'No files for subject %d in baseFolder. Skipping.', sid);
                        continue;
                    end
                end

                % After a bootstrap every indexed file was just accumulated —
                % nothing new to check, skip the candidate loop entirely.
                if wasBootstrapped
                    fprintf('  All files accumulated during bootstrap -- nothing to patch.\n');
                    continue;
                end

                % ── Gather candidate individual combo files ─────────────
                if ~obj.subjectMap.isKey(sid)
                    fprintf('  No new files indexed for subject %d -- skipping.\n', sid);
                    continue;
                end
                candidateFiles = obj.subjectMap(sid);

                % ── Check raw file present in output folder ─────────────
                accRawFname = sprintf('S%d_%s_raw.mat', sid, obj.algorithmName);
                accRawFpath = fullfile(obj.outputFolder, accRawFname);
                if ~isfile(accRawFpath)
                    srcRawFpath = fullfile(obj.baseFolder, accRawFname);
                    if isfile(srcRawFpath)
                        loaded2 = load(srcRawFpath, 'raw');
                        raw = struct( ...
                            'calibration', loaded2.raw.cal,    ...
                            'closed',      loaded2.raw.closed, ...
                            'open',        loaded2.raw.open);  %#ok<NASGU>
                        save(accRawFpath, 'raw', '-v7.3');
                        fprintf('  Copied raw signals -> %s\n', accRawFname);
                    else
                        warning('accumulateResults:patch', ...
                            'Raw file not found for S%d -- raw signals will be missing.', sid);
                    end
                end

                % ── Build a fast param-key set from existing meta ────────
                % Keys are built from param structs without loading any signal data.
                existingKeys = containers.Map('KeyType','char','ValueType','logical');
                for ki = 1:existingN
                    if ~isempty(existingParams{ki})
                        existingKeys(obj.paramStructToKey(existingParams{ki})) = true;
                    end
                end

                % ── Identify and append only new combos ─────────────────
                % Duplicate check uses only the config struct (tiny), not the
                % signal data. parseFile() is only called for genuinely new files.
                nextK    = existingN + 1;
                newAdded = 0;

                for fi = 1:numel(candidateFiles)
                    fpath = candidateFiles{fi};

                    % Fast check: load only the config from the file
                    try
                        S = load(fpath, 'results');
                        candidateParams = S.results.config.algorithm.parameters;
                    catch e
                        warning('accumulateResults:patch', ...
                            'Cannot read config from %s: %s', fpath, e.message);
                        continue;
                    end

                    candidateKey = obj.paramStructToKey(candidateParams);
                    if existingKeys.isKey(candidateKey)
                        continue;   % already accumulated — skip without loading signals
                    end

                    % Genuinely new — now do the full parse and write
                    try
                        entry = obj.parseFile(fpath);
                    catch e
                        warning('accumulateResults:patch', ...
                            'Skipping %s: %s', fpath, e.message);
                        continue;
                    end

                    comboFields = {'parameters','timestamp', ...
                                   'cleanCalibration','cleanClosed','cleanOpen', ...
                                   'modifiedMask','timeProbe','probeRaw','probeClean'};
                    combo = struct();                               %#ok<NASGU>
                    for f = 1:numel(comboFields)
                        combo.(comboFields{f}) = entry.(comboFields{f});
                    end
                    comboFname = sprintf('S%d_%s_combo_%d.mat', sid, obj.algorithmName, nextK);
                    comboFpath = fullfile(obj.outputFolder, comboFname);
                    save(comboFpath, 'combo', '-v7.3');

                    existingParams{nextK}    = entry.parameters;
                    existingKeys(candidateKey) = true;
                    fprintf('  [+] combo %d  <- %s\n', nextK, obj.paramStructToStr(entry.parameters));

                    nextK    = nextK    + 1;
                    newAdded = newAdded + 1;
                end

                % ── Rewrite _meta.mat if anything was added ───────────────
                if newAdded > 0
                    nCombos    = nextK - 1;                         %#ok<NASGU>
                    paramsList = existingParams;                     %#ok<NASGU>
                    save(metaFpath, 'nCombos', 'paramsList', '-v7.3');
                    fprintf('  meta updated: %d total combos (%d new).\n', nCombos, newAdded);
                else
                    fprintf('  No new combos found.\n');
                end

                nAdded = nAdded + newAdded;
            end

            fprintf('\npatch done. %d new combos added across %d subject(s).\n', ...
                nAdded, numel(targetIDs));
        end

        % ================================================================
        % patchCleanup — delete individual run files for already-patched combos
        % ================================================================
        function patchCleanup(obj, subjectID)
            % PATCHCLEANUP - Delete individual combo run files that are already
            %   present in the split accumulated layout.
            %   Safe: only deletes a file if its params are found in _meta.mat.
            %
            %   patchCleanup()            -- all subjects
            %   patchCleanup(subjectID)   -- one subject only
            %
            %   Example:
            %     acc.patch();
            %     acc.patchCleanup();   % clean up only the ingested files

            if nargin < 2
                targetIDs = obj.subjects;
            else
                targetIDs = subjectID(:)';
            end

            % Auto-index if not yet done
            if obj.subjectMap.Count == 0
                obj.locateFiles();
            end

            nDeleted = 0;

            for s = 1:numel(targetIDs)
                sid = targetIDs(s);

                metaFname = sprintf('S%d_%s_meta.mat', sid, obj.algorithmName);
                metaFpath = fullfile(obj.outputFolder, metaFname);
                if ~isfile(metaFpath)
                    fprintf('  S%d: no meta file -- skipping cleanup.\n', sid);
                    continue;
                end

                loaded       = load(metaFpath, 'paramsList');
                knownParams  = loaded.paramsList;

                % Build fast key set
                knownKeys = containers.Map('KeyType','char','ValueType','logical');
                for ki = 1:numel(knownParams)
                    if ~isempty(knownParams{ki})
                        knownKeys(obj.paramStructToKey(knownParams{ki})) = true;
                    end
                end

                if ~obj.subjectMap.isKey(sid), continue; end
                filePaths = obj.subjectMap(sid);

                for fi = 1:numel(filePaths)
                    fpath = filePaths{fi};
                    try
                        S = load(fpath, 'results');
                        candidateKey = obj.paramStructToKey(S.results.config.algorithm.parameters);
                    catch
                        continue;
                    end
                    if knownKeys.isKey(candidateKey)
                        delete(fpath);
                        nDeleted = nDeleted + 1;
                    end
                end
            end

            fprintf('patchCleanup done. %d files deleted.\n', nDeleted);
        end

        % ================================================================
        % listCombos — print a table of all accumulated combos
        % ================================================================
        function listCombos(obj, subjectID)
            % LISTCOMBOS - Print all parameter combos in the split meta file.
            %   Useful for checking what's already been accumulated.
            %
            %   listCombos()            -- all subjects
            %   listCombos(subjectID)   -- one subject only

            if nargin < 2
                targetIDs = obj.subjects;
            else
                targetIDs = subjectID(:)';
            end

            for s = 1:numel(targetIDs)
                sid = targetIDs(s);
                metaFname = sprintf('S%d_%s_meta.mat', sid, obj.algorithmName);
                metaFpath = fullfile(obj.outputFolder, metaFname);
                if ~isfile(metaFpath)
                    fprintf('S%d: no meta file found.\n', sid);
                    continue;
                end
                loaded = load(metaFpath, 'nCombos', 'paramsList');
                fprintf('\nS%d  %s  --  %d combos:\n', sid, obj.algorithmName, loaded.nCombos);
                for k = 1:loaded.nCombos
                    fprintf('  [%3d]  %s\n', k, obj.paramStructToStr(loaded.paramsList{k}));
                end
            end
        end

        % ================================================================
        % resplit — migrate monolithic file(s) to split layout
        % ================================================================
        function resplit(obj, subjectID)
            % RESPLIT - Convert monolithic _accumulated.mat to split layout.
            %   Reads the existing monolithic file and writes the three split
            %   files (_raw.mat, _meta.mat, _combo_k.mat × N). The monolithic
            %   file is left intact (delete manually once verified).
            %
            %   resplit()            -- migrate all subjects in obj.subjects
            %   resplit(subjectID)   -- migrate one subject only
            %
            %   Example:
            %     acc = accumulateResults(folder, 'vanilla-asr', [1 2 3]);
            %     acc.resplit();        % migrate all three subjects

            if nargin < 2
                targetIDs = obj.subjects;
            else
                targetIDs = subjectID;
            end

            for s = 1:numel(targetIDs)
                sid = targetIDs(s);

                mono = sprintf('S%d_%s_accumulated.mat', sid, obj.algorithmName);
                mpath = fullfile(obj.outputFolder, mono);
                if ~isfile(mpath)
                    warning('accumulateResults:resplit', ...
                        'Monolithic file not found for subject %d -- skipping.\n  Expected: %s', ...
                        sid, mpath);
                    continue;
                end

                fprintf('  resplit S%d %s...\n', sid, obj.algorithmName);
                loaded          = load(mpath, 'raw', 'subject_results');
                raw_loaded      = loaded.raw;
                sr_all          = loaded.subject_results;
                nC              = numel(sr_all);

                obj.saveSplit(sid, raw_loaded, sr_all, nC);
                fprintf('    -> %d split files written.\n', nC + 2);
            end
        end

        % ================================================================
        % Delete individual run files after accumulation
        % ================================================================
        function cleanup(obj, subjectID)
            % CLEANUP - Delete individual per-run .mat files from baseFolder.
            %   Checks that an accumulated output exists before deleting anything.
            %   Accepts both monolithic (_accumulated.mat) and split (_meta.mat) layouts.
            %   Calls locateFiles() automatically if subjectMap is not yet populated.
            %
            %   cleanup()          -- delete files for ALL subjects
            %   cleanup(subjectID) -- delete files for one subject only

            if nargin < 2
                targetIDs = obj.subjects;
            else
                targetIDs = subjectID;
            end

            % Auto-index if not yet done (e.g. calling cleanup standalone)
            if obj.subjectMap.Count == 0
                obj.locateFiles();
            end

            nDeleted = 0;

            for s = 1:numel(targetIDs)
                sid = targetIDs(s);

                % Safety: confirm accumulated output exists (monolithic OR split) before deleting
                monoFname = sprintf('S%d_%s_accumulated.mat', sid, obj.algorithmName);
                metaFname = sprintf('S%d_%s_meta.mat',        sid, obj.algorithmName);
                monoOk    = isfile(fullfile(obj.outputFolder, monoFname));
                metaOk    = isfile(fullfile(obj.outputFolder, metaFname));
                if ~monoOk && ~metaOk
                    warning('accumulateResults:cleanup', ...
                        'No accumulated output for subject %d -- skipping cleanup.\n  Checked: %s  /  %s', ...
                        sid, monoFname, metaFname);
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

        function saveSplit(obj, subjectID, raw, subject_results, nValid)
            % SAVESPLIT - Write the three split-layout files for one subject.
            %   _raw.mat      : raw struct  (.calibration/.closed/.open)
            %   _meta.mat     : nCombos + parameters cell array (no signal data)
            %   _combo_k.mat  : one file per combo with cleaned signals + diagnostics

            % ── _raw.mat ──────────────────────────────────────────────
            rawFname = sprintf('S%d_%s_raw.mat', subjectID, obj.algorithmName);
            rawFpath = fullfile(obj.outputFolder, rawFname);
            save(rawFpath, 'raw', '-v7.3');

            % ── _meta.mat ─────────────────────────────────────────────
            nCombos = nValid;                                      %#ok<NASGU>
            paramsList = cell(1, nValid);
            for k = 1:nValid
                paramsList{k} = subject_results(k).parameters;
            end
            metaFname = sprintf('S%d_%s_meta.mat', subjectID, obj.algorithmName);
            metaFpath = fullfile(obj.outputFolder, metaFname);
            save(metaFpath, 'nCombos', 'paramsList', '-v7.3');

            % ── _combo_k.mat (one per combo) ──────────────────────────
            comboFields = {'parameters','timestamp', ...
                           'cleanCalibration','cleanClosed','cleanOpen', ...
                           'modifiedMask','timeProbe','probeRaw','probeClean'};
            for k = 1:nValid
                combo = struct();                                   %#ok<NASGU>
                for f = 1:numel(comboFields)
                    combo.(comboFields{f}) = subject_results(k).(comboFields{f});
                end
                comboFname = sprintf('S%d_%s_combo_%d.mat', subjectID, obj.algorithmName, k);
                comboFpath = fullfile(obj.outputFolder, comboFname);
                save(comboFpath, 'combo', '-v7.3');
            end

            fprintf('    Saved (split): %s  [%d combos]\n', ...
                sprintf('S%d_%s_{raw,meta,combo_k}', subjectID, obj.algorithmName), nValid);
        end

        function found = paramExistsInList(~, params, paramsList)
            % PARAMEXISTSINLIST - Returns true if params matches any entry in paramsList.
            %   Comparison is field-by-field, numeric values compared with isequal.
            %   String fields compared case-insensitively.
            found = false;
            fields = fieldnames(params);
            for k = 1:numel(paramsList)
                if isempty(paramsList{k}), continue; end
                candidate = paramsList{k};
                match = true;
                for f = 1:numel(fields)
                    fn = fields{f};
                    if ~isfield(candidate, fn)
                        match = false; break;
                    end
                    a = params.(fn);
                    b = candidate.(fn);
                    if ischar(a) || isstring(a)
                        if ~strcmpi(char(a), char(b)), match = false; break; end
                    else
                        if ~isequal(a, b), match = false; break; end
                    end
                end
                if match, found = true; return; end
            end
        end

        function key = paramStructToKey(~, params)
            % PARAMSTRUCTTOKEY - Deterministic string key from a param struct.
            %   Used for O(1) duplicate detection via containers.Map.
            %   Fields are sorted alphabetically so key is order-independent.
            fields = sort(fieldnames(params));
            parts  = cell(1, numel(fields));
            for f = 1:numel(fields)
                v = params.(fields{f});
                if ischar(v) || isstring(v)
                    parts{f} = sprintf('%s=%s', fields{f}, char(v));
                elseif isscalar(v)
                    parts{f} = sprintf('%s=%.10g', fields{f}, v);
                else
                    parts{f} = sprintf('%s=[%s]', fields{f}, num2str(v(:)', '%.10g '));
                end
            end
            key = strjoin(parts, '|');
        end

        function str = paramStructToStr(~, params)
            % PARAMSTRUCTTOSTR - One-line summary of a parameter struct.
            fields = fieldnames(params);
            parts  = cell(1, numel(fields));
            for f = 1:numel(fields)
                v = params.(fields{f});
                if ischar(v) || isstring(v)
                    parts{f} = sprintf('%s=%s', fields{f}, char(v));
                elseif isscalar(v)
                    parts{f} = sprintf('%s=%g', fields{f}, v);
                else
                    parts{f} = sprintf('%s=[%s]', fields{f}, num2str(v));
                end
            end
            str = strjoin(parts, '  ');
        end

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
