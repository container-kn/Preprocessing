%% run_analysis_batch.m
%
% Batch analysis script — scans a directory for accumulated .mat files,
% runs every analysis class on each (subject × algorithm) pair, and saves
% results + summary tables.
%
% ── WHAT IT DOES ──────────────────────────────────────────────────────────
%   For every S<id>_<algorithm>_accumulated.mat found in accumulatedDir it:
%
%   1.  Loads the file via ExperimentAnalysis
%   2.  Runs all enabled analyses:
%         a. BlinkAnalysis      — blink detection and reduction metrics
%         b. SignalStatistics   — RMS, SNR, RRMSE, correlation per segment
%         c. QuietRegionAnalysis— metrics on blink-free regions (needs blink)
%         d. SpectralAnalysis   — band power before/after per segment
%         e. TemporalAnalysis   — windowed signal + subspace + compaction + drift
%   3.  Saves a per-subject .mat with all analysis objects
%   4.  Appends one summary row per combo to a cross-subject table
%   5.  Writes the cross-subject table as a .csv and .mat at the end
%
% ── HOW TO RUN ────────────────────────────────────────────────────────────
%   1. Edit the CONFIGURATION section below
%   2. Run the whole script:  run_analysis_batch
%   3. Results land in outputDir
%
% ── OUTPUT FILES ──────────────────────────────────────────────────────────
%   <outputDir>/
%     analysis_S<id>_<algo>.mat   — all analysis objects for one subject
%     summary_table.csv           — all subjects × combos × metrics (flat)
%     summary_table.mat           — same as MATLAB table
%     batch_log.txt               — per-run timing and error log
%
% ── DEPENDENCIES ──────────────────────────────────────────────────────────
%   Requires all analysis class files on the path:
%     ExperimentAnalysis, BlinkAnalysis, SignalStatistics,
%     QuietRegionAnalysis, SpectralAnalysis, TemporalAnalysis
%   Optional (graceful skip if absent):
%     ICAAnalysis  — requires EEGLAB (runica + iclabel)
%     probeAnalysis
%   EEGLAB (for SpectralAnalysis welch, ICA, robust cov):  optional
%
% =========================================================================

%% ── CONFIGURATION ────────────────────────────────────────────────────────
% Edit this section only.

% --- Paths ---
accumulatedDir = 'K:\AccumulatedResults\GraphASR_Results';    % folder with S*_*_accumulated.mat
outputDir      = 'K:\Analysis\GraphASR';   % where results will be written

% --- Sampling rate ---
fs = 500;   % Hz — must match recordings

% --- Which analyses to run ---
RUN_BLINKS    = true;
RUN_SIGNAL    = true;
RUN_QUIET     = true;   % needs RUN_BLINKS = true
RUN_SPECTRAL  = true;
RUN_TEMPORAL  = true;
RUN_ICA       = false;  % expensive — requires EEGLAB; set true to enable

% --- BlinkAnalysis options ---
blinkOpts = { ...
    'method',     'mad',  ... % 'mad' | 'mean'
    'multiplier',  5,     ...
    'bandpass',   [1 10]};

% --- SpectralAnalysis: no options (compute() takes none) ---

% --- TemporalAnalysis options ---
temporalOpts = { ...
    'windowSec',         2,    ...
    'overlap',           0.5,  ...  % 50% overlap (alternative to hopSec)
    'subspaceRank',      10,   ...
    'compactionThresh',  0.9,  ...
    'kRef',              [],   ...  % [] = same as subspaceRank
    'useRobustCov',      true};

% --- ICAAnalysis options (only used when RUN_ICA=true) ---
icaOpts = { ...
    'segment',            'closed+open', ...
    'maxSamples',         30000,         ...
    'preservationThresh', 0.8};

% --- Summary metrics to extract for the cross-subject table ---
% Format: {analyser, fieldname, segment}
%   analyser  : 'blink' | 'signalStats' | 'temporal' | 'spectral'
%   fieldname : passed to that analyser's getSummary()
%   segment   : 'closed' | 'open' (ignored for blink)
SUMMARY_METRICS = { ...
    % Blink
    'blink', 'blinkReductionPct_closed',    '';
    'blink', 'blinkReductionPct_open',       '';
    'blink', 'persistenceRate_closed',       '';
    'blink', 'persistenceRate_open',         '';
    'blink', 'removalRate_closed',           '';
    % Signal
    'signalStats', 'rrmse',        'closed';
    'signalStats', 'corr',         'closed';
    'signalStats', 'snrDB',        'closed';
    'signalStats', 'signalChange', 'closed';
    'signalStats', 'rrmse',        'open';
    'signalStats', 'corr',         'open';
    % Temporal
    'temporal', 'signal.rrmse_mean',                    'closed';
    'temporal', 'signal.corr_mean',                     'closed';
    'temporal', 'signal.rrmse_p95',                     'closed';
    'temporal', 'subspace.angle_mean_mean',              'closed';
    'temporal', 'subspace.energy_removed_mean',          'closed';
    'temporal', 'subspace.cov_fro_norm_mean',            'closed';
    'temporal', 'energyCompaction.dimReduction_mean',    'closed';
    'temporal', 'energyCompaction.kCompact_raw_mean',    'closed';
    'temporal', 'energyCompaction.kCompact_clean_mean',  'closed';
    'temporal', 'calibrationDrift.energyCaptureGain_mean','closed';
    'temporal', 'calibrationDrift.angleFromRefGain_mean', 'closed';
    % Spectral
    'spectral', 'totalPowerChange',  'closed';
    'spectral', 'spectralCorr',      'closed';
    'spectral', 'alpha.totalPowerAfter', 'closed';
    'spectral', 'beta.totalPowerAfter',  'closed';
    'spectral', 'gamma.totalPowerAfter', 'closed';
};

% =========================================================================
%% ── SETUP ────────────────────────────────────────────────────────────────

if ~isfolder(outputDir), mkdir(outputDir); end

logFile  = fullfile(outputDir, 'batch_log.txt');
logFid   = fopen(logFile, 'w');
logPrint(logFid, '=== run_analysis_batch ===');
logPrint(logFid, 'Started: %s', datetime('now'));
logPrint(logFid, 'Input  : %s', accumulatedDir);
logPrint(logFid, 'Output : %s', outputDir);

% --- Discover accumulated files ---
d = dir(fullfile(accumulatedDir, 'S*_*_accumulated.mat'));

if isempty(d)
    error('run_analysis_batch:noFiles', ...
        'No S*_*_accumulated.mat files found in:\n  %s', accumulatedDir);
end

logPrint(logFid, 'Found %d accumulated files\n', numel(d));
fprintf('\n=== run_analysis_batch ===\n');
fprintf('Found %d accumulated file(s) in:\n  %s\n\n', numel(d), accumulatedDir);

% --- Master summary table accumulator ---
allRows = {};   % cell array of tables, one per (subject x algo) pair

tBatch = tic;

%% ── MAIN LOOP ────────────────────────────────────────────────────────────

for fi = 1:numel(d)

    fname = d(fi).name;
    fpath = fullfile(d(fi).folder, fname);

    % ── Parse filename: S<id>_<algo>_accumulated.mat ─────────────────
    tok = regexp(fname, '^S(\d+)_(.+)_accumulated\.mat$', 'tokens');
    if isempty(tok)
        logPrint(logFid, '[SKIP] Cannot parse filename: %s', fname);
        continue;
    end
    subjectID  = str2double(tok{1}{1});
    algoName   = tok{1}{2};

    fprintf('[%d/%d] S%d  %s\n', fi, numel(d), subjectID, algoName);
    logPrint(logFid, '\n[%d/%d] S%d  %s', fi, numel(d), subjectID, algoName);

    tSubject = tic;

    % ── Load via ExperimentAnalysis ───────────────────────────────────
    try
        ana = ExperimentAnalysis(accumulatedDir, algoName, subjectID, fs);
        ana.load();
    catch e
        logPrint(logFid, '  ERROR loading: %s', e.message);
        fprintf('  ERROR loading S%d %s: %s\n', subjectID, algoName, e.message);
        continue;
    end

    fprintf('  Loaded: %d combo(s)\n', ana.nCombos);

    % ── a. BlinkAnalysis ─────────────────────────────────────────────
    if RUN_BLINKS
        try
            t0 = tic;
            ana.computeBlinks(blinkOpts{:});
            logPrint(logFid, '  blink       %.1fs', toc(t0));
            fprintf('  [ok] blink (%.1fs)\n', toc(t0));
        catch e
            logPrint(logFid, '  blink  ERROR: %s', e.message);
            fprintf('  [fail] blink: %s\n', e.message);
        end
    end

    % ── b. SignalStatistics ───────────────────────────────────────────
    if RUN_SIGNAL
        try
            t0 = tic;
            ana.computeSignalStats();
            logPrint(logFid, '  signal      %.1fs', toc(t0));
            fprintf('  [ok] signal (%.1fs)\n', toc(t0));
        catch e
            logPrint(logFid, '  signal ERROR: %s', e.message);
            fprintf('  [fail] signal: %s\n', e.message);
        end
    end

    % ── c. QuietRegionAnalysis ────────────────────────────────────────
    if RUN_QUIET && RUN_BLINKS && ~isempty(ana.blink)
        try
            t0 = tic;
            qr = QuietRegionAnalysis( ...
                ana.raw, ana.subject_results, ana.fs, ana.blink, ...
                ana.subjectID, ana.algorithmName);
            qr.compute();
            ana.quietRegion = qr;   %#ok<STRNU> — store on ana for saving
            logPrint(logFid, '  quiet       %.1fs', toc(t0));
            fprintf('  [ok] quiet region (%.1fs)\n', toc(t0));
        catch e
            logPrint(logFid, '  quiet  ERROR: %s', e.message);
            fprintf('  [fail] quiet region: %s\n', e.message);
        end
    elseif RUN_QUIET && (~RUN_BLINKS || isempty(ana.blink))
        logPrint(logFid, '  quiet  SKIP (need blink first)');
        fprintf('  [skip] quiet region (need blink)\n');
    end

    % ── d. SpectralAnalysis ───────────────────────────────────────────
    if RUN_SPECTRAL
        try
            t0 = tic;
            sp = SpectralAnalysis( ...
                ana.raw, ana.subject_results, ana.fs, ...
                ana.subjectID, ana.algorithmName);
            sp.compute();
            ana.spectral = sp;   %#ok<STRNU>
            logPrint(logFid, '  spectral    %.1fs', toc(t0));
            fprintf('  [ok] spectral (%.1fs)\n', toc(t0));
        catch e
            logPrint(logFid, '  spectral ERROR: %s', e.message);
            fprintf('  [fail] spectral: %s\n', e.message);
        end
    end

    % ── e. TemporalAnalysis ───────────────────────────────────────────
    if RUN_TEMPORAL
        try
            t0 = tic;
            ta = TemporalAnalysis( ...
                ana.raw, ana.subject_results, ana.fs, ...
                ana.subjectID, ana.algorithmName);
            ta.compute(temporalOpts{:});
            ana.temporal = ta;   %#ok<STRNU>
            logPrint(logFid, '  temporal    %.1fs', toc(t0));
            fprintf('  [ok] temporal (%.1fs)\n', toc(t0));
        catch e
            logPrint(logFid, '  temporal ERROR: %s', e.message);
            fprintf('  [fail] temporal: %s\n', e.message);
        end
    end

    % ── f. ICAAnalysis (optional, expensive) ─────────────────────────
    if RUN_ICA
        try
            t0 = tic;
            ana.computeICA(icaOpts{:});
            logPrint(logFid, '  ica         %.1fs', toc(t0));
            fprintf('  [ok] ICA (%.1fs)\n', toc(t0));
        catch e
            logPrint(logFid, '  ica    ERROR: %s', e.message);
            fprintf('  [fail] ICA: %s\n', e.message);
        end
    end

    % ── Save per-subject analysis .mat ───────────────────────────────
    saveName = sprintf('analysis_S%d_%s.mat', subjectID, algoName);
    savePath = fullfile(outputDir, saveName);
    try
        save(savePath, 'ana', '-v7.3');
        logPrint(logFid, '  saved: %s', saveName);
        fprintf('  [saved] %s\n', saveName);
    catch e
        logPrint(logFid, '  SAVE ERROR: %s', e.message);
        fprintf('  [fail] save: %s\n', e.message);
    end

    % ── Build summary rows for cross-subject table ────────────────────
    rows = buildSummaryRows(ana, subjectID, algoName, SUMMARY_METRICS);
    if ~isempty(rows)
        allRows{end+1} = rows; %#ok<AGROW>
    end

    elapsed = toc(tSubject);
    logPrint(logFid, '  Total: %.1fs', elapsed);
    fprintf('  Done (%.1fs)\n\n', elapsed);

end % file loop

%% ── WRITE CROSS-SUBJECT TABLE ────────────────────────────────────────────

if ~isempty(allRows)
    try
        summaryTable = vertcat(allRows{:});

        % CSV
        csvPath = fullfile(outputDir, 'summary_table.csv');
        writetable(summaryTable, csvPath);

        % MAT
        matPath = fullfile(outputDir, 'summary_table.mat');
        save(matPath, 'summaryTable');

        fprintf('Summary table: %d rows x %d columns\n', ...
            height(summaryTable), width(summaryTable));
        fprintf('  -> %s\n', csvPath);

        logPrint(logFid, '\nSummary: %d rows x %d cols -> %s', ...
            height(summaryTable), width(summaryTable), csvPath);
    catch e
        fprintf('[fail] Summary table: %s\n', e.message);
        logPrint(logFid, 'Summary table ERROR: %s', e.message);
    end
else
    fprintf('No summary rows collected.\n');
end

totalElapsed = toc(tBatch);
logPrint(logFid, '\n=== Batch complete: %.1fs ===', totalElapsed);
fprintf('\n=== Batch complete in %.1fs ===\n', totalElapsed);
fprintf('Output: %s\n\n', outputDir);
fclose(logFid);


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

% ── buildSummaryRows ─────────────────────────────────────────────────────
function rows = buildSummaryRows(ana, subjectID, algoName, metrics)
% BUILDSUMMARYROWS
%   Returns a table with one row per parameter combo.
%   Columns: subjectID, algorithmName, comboIdx, all scalar parameters,
%   then one column per entry in the metrics cell array.

    rows = table();

    if ana.nCombos == 0
        return;
    end

    % --- Base columns: identity ---
    nC = ana.nCombos;
    rows.subjectID    = repmat(subjectID, nC, 1);
    rows.algorithmName = repmat({algoName}, nC, 1);
    rows.comboIdx     = (1:nC)';

    % --- Parameter columns (from first combo, scalar fields) ---
    if ~isempty(ana.subject_results) && isfield(ana.subject_results(1), 'parameters')
        params = ana.subject_results(1).parameters;
        pf = fieldnames(params);
        for pi = 1:numel(pf)
            colName = pf{pi};
            vals = NaN(nC, 1);
            for c = 1:nC
                v = ana.subject_results(c).parameters.(colName);
                if isscalar(v) && isnumeric(v)
                    vals(c) = v;
                end
            end
            rows.(colName) = vals;
        end
    end

    % --- Metric columns ---
    for mi = 1:size(metrics, 1)
        analyser  = metrics{mi, 1};
        fieldName = metrics{mi, 2};
        segment   = metrics{mi, 3};

        % Safe column name (replace . and spaces)
        colName = [analyser, '_', strrep(strrep(fieldName, '.', '_'), ' ', '_')];
        if ~isempty(segment)
            colName = [colName, '_', segment]; %#ok<AGROW>
        end

        vals = NaN(nC, 1);
        try
            switch analyser
                case 'blink'
                    if ~isempty(ana.blink)
                        vals = ana.blink.getSummary(fieldName)';
                    end
                case 'signalStats'
                    if ~isempty(ana.signalStats)
                        % getSummary(field, segment, signalType)
                        % default signalType = 'diff'
                        v = ana.signalStats.getSummary(fieldName, segment, 'diff');
                        vals = v(:);
                    end
                case 'spectral'
                    if isfield(ana, 'spectral') && ~isempty(ana.spectral)
                        v = ana.spectral.getSummary(fieldName, segment);
                        vals = v(:);
                    end
                case 'temporal'
                    if isfield(ana, 'temporal') && ~isempty(ana.temporal)
                        v = ana.temporal.getSummary(fieldName, segment);
                        vals = v(:);
                    end
                case 'ica'
                    if ~isempty(ana.ica)
                        v = ana.ica.getSummary(fieldName);
                        vals = v(:);
                    end
            end
        catch
            % leave NaN — don't crash on missing field
        end

        rows.(colName) = vals;
    end
end

% ── logPrint ──────────────────────────────────────────────────────────────
function logPrint(fid, fmt, varargin)
% LOGPRINT - Print to log file with timestamp.
    try
        msg = sprintf(fmt, varargin{:});
        fprintf(fid, '[%s] %s\n', datestr(now, 'HH:MM:SS'), msg);
    catch
    end
end
