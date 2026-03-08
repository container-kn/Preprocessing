%% run_analysis_batch.m
%
% Batch analysis script — streams one combo at a time, analyses it, saves
% results incrementally. Never loads more than batchSize combos into memory.
%
% ── STREAMING WORKFLOW ────────────────────────────────────────────────────
%   For each (subject × algorithm) pair:
%     1. loadMeta()  — load raw signals + combo count only
%     2. For each batch of combos:
%          loadCombo(k)       — load cleaned signals for combo k only
%          run all analyses   — on this batch only
%          save/merge results — append to per-subject .mat
%          unloadCombo(k)     — free signal memory before next batch
%     3. Write summary table row after all combos done
%
% ── OUTPUT FILES ──────────────────────────────────────────────────────────
%   <outputDir>/
%     analysis_S<id>_<algo>.mat   — comboResults struct + doneFlags
%     summary_table.csv / .mat    — all subjects × combos × metrics
%     batch_log.txt               — timing and error log (append mode)
%
% ── RESUMING ──────────────────────────────────────────────────────────────
%   RESUME = true  skips combos already in analysis_S<id>_<algo>.mat
%   Safe to re-run after a crash or interruption.
%
% =========================================================================

%% ── CONFIGURATION ────────────────────────────────────────────────────────

accumulatedDir = 'J:\ASRProjct\results\accumulated';
outputDir      = 'J:\ASRProjct\analysis\graph-asr';
fs             = 500;

% Number of combos to process per batch before saving and freeing memory.
% 1 = minimum RAM.  5-10 = good balance.  Inf = load all at once.
batchSize = 5;

% Skip combos whose results are already saved (safe after crash)
RESUME = true;

% Which analyses to run
RUN_BLINKS   = true;
RUN_SIGNAL   = true;
RUN_QUIET    = true;    % requires RUN_BLINKS
RUN_SPECTRAL = true;
RUN_TEMPORAL = true;
RUN_ICA      = false;

% Analysis options
blinkOpts    = {'method','mad','multiplier',5,'bandpass',[1 10]};
spectralOpts = {};
quietOpts    = {};
temporalOpts = {'windowSec',2,'overlap',0.5,'subspaceRank',10, ...
                'compactionThresh',0.9,'useRobustCov',true};
icaOpts      = {'segment','closed+open','maxSamples',30000,'preservationThresh',0.8};

% Summary metrics: {metricName, fieldName, segment}
% signalStats fieldName format: '<segment>.<signalType>.<field>'
%   e.g. 'closed.diff.rrmse'  'closed.diff.corr'  'open.diff.rrmse'
% blink:       'blinkReductionPct_closed'  'persistenceRate_closed' etc.
% spectral:    'alpha.totalPowerAfter'  'totalPowerChange' etc.
% temporal:    'signal.rrmse_mean'  'subspace.angle_mean_mean' etc.
% quietRegion: 'corr_median'  'rre_median'  'quietFraction' etc.
SUMMARY_METRICS = { ...
    'blink',        'blinkReductionPct_closed',              '';
    'blink',        'blinkReductionPct_open',                '';
    'blink',        'persistenceRate_closed',                '';
    'blink',        'removalRate_closed',                    '';
    'signalStats',  'closed.diff.rrmse',                     '';
    'signalStats',  'closed.diff.corr',                      '';
    'signalStats',  'open.diff.rrmse',                       '';
    'quietRegion',  'corr_median',                           'closed';
    'quietRegion',  'rre_median',                            'closed';
    'quietRegion',  'quietFraction',                         'closed';
    'spectral',     'totalPowerChange',                      'closed';
    'spectral',     'alpha.totalPowerAfter',                 'closed';
    'spectral',     'alpha.relChange',                       'closed';
    'temporal',     'signal.rrmse_mean',                     'closed';
    'temporal',     'signal.corr_mean',                      'closed';
    'temporal',     'subspace.angle_mean_mean',              'closed';
    'temporal',     'energyCompaction.dimReduction_mean',    'closed';
    'temporal',     'calibrationDrift.energyCaptureGain_mean','closed';
};

% =========================================================================
%% ── SETUP ────────────────────────────────────────────────────────────────

if ~isfolder(outputDir), mkdir(outputDir); end

logFile = fullfile(outputDir, 'batch_log.txt');
logFid  = fopen(logFile, 'a');
logPrint(logFid, '=== run_analysis_batch (streaming) ===');
logPrint(logFid, 'Started: %s | batchSize=%d | resume=%d', datetime('now'), batchSize, RESUME);

dMeta   = dir(fullfile(accumulatedDir, 'S*_*_meta.mat'));
dMono   = dir(fullfile(accumulatedDir, 'S*_*_accumulated.mat'));
entries = discoverSubjects(dMeta, '_meta.mat', dMono, '_accumulated.mat');

if isempty(entries)
    error('run_analysis_batch:noFiles', 'No accumulated files found in:\n  %s', accumulatedDir);
end

fprintf('\n=== run_analysis_batch ===\n');
fprintf('Found %d pair(s)  |  batchSize=%d  |  resume=%d\n\n', numel(entries), batchSize, RESUME);

allSummaryRows = {};
tBatch = tic;

%% ── MAIN LOOP ────────────────────────────────────────────────────────────

for fi = 1:numel(entries)

    subjectID = entries(fi).subjectID;
    algoName  = entries(fi).algoName;
    isSplit   = entries(fi).isSplit;

    fprintf('[%d/%d] S%d  %s\n', fi, numel(entries), subjectID, algoName);
    logPrint(logFid, '\n[%d/%d] S%d  %s', fi, numel(entries), subjectID, algoName);

    saveFile = fullfile(outputDir, sprintf('analysis_S%d_%s.mat', subjectID, algoName));
    tSubject = tic;

    % ── Load raw + meta (no signal data yet) ──────────────────────────
    try
        ana = ExperimentAnalysis(accumulatedDir, algoName, subjectID, fs);
        if isSplit
            ana.loadMeta();
        else
            ana.load();   % monolithic: loads all at once
        end
    catch e
        logPrint(logFid, '  ERROR loading: %s', e.message);
        fprintf('  ERROR: %s\n', e.message);
        continue;
    end
    fprintf('  %d combos total\n', ana.nCombos);

    % ── Load previous partial results if resuming ─────────────────────
    comboResults = initComboResults(ana.nCombos);
    doneFlags    = false(1, ana.nCombos);

    if RESUME && isfile(saveFile)
        try
            saved        = load(saveFile, 'comboResults', 'doneFlags');
            comboResults = saved.comboResults;
            doneFlags    = saved.doneFlags;
            fprintf('  Resuming: %d/%d already done\n', sum(doneFlags), ana.nCombos);
            logPrint(logFid, '  Resuming: %d/%d done', sum(doneFlags), ana.nCombos);
        catch
            fprintf('  No previous results found -- starting fresh.\n');
        end
    end

    combosToRun = find(~doneFlags);
    if isempty(combosToRun)
        fprintf('  All combos already done -- skipping.\n');
        rows = buildSummaryRows(comboResults, subjectID, algoName, ana.nCombos, SUMMARY_METRICS);
        if ~isempty(rows), allSummaryRows{end+1} = rows; end %#ok<AGROW>
        continue;
    end

    % ── Stream through combos in batches ──────────────────────────────
    for bStart = 1 : batchSize : numel(combosToRun)

        bEnd    = min(bStart + batchSize - 1, numel(combosToRun));
        batch_k = combosToRun(bStart:bEnd);

        fprintf('  Batch combos [%s] (%d of %d remaining)...\n', ...
            num2str(batch_k), numel(combosToRun) - bStart + 1, numel(combosToRun));

        % Load signal data for this batch
        loadedOk = true(size(batch_k));   % track which combos loaded cleanly
        if isSplit
            for bi = 1:numel(batch_k)
                k = batch_k(bi);
                try
                    ana.loadCombo(k);
                catch e
                    fprintf('  [fail] loadCombo(%d): %s\n', k, e.message);
                    loadedOk(bi) = false;
                end
            end
        end

        % Skip the entire batch if nothing loaded
        if ~any(loadedOk)
            fprintf('  [skip] entire batch — all loads failed\n');
            continue;
        end

        % Restrict to successfully-loaded combos only
        batch_k_ok = batch_k(loadedOk);

        % Build a batch-restricted ExperimentAnalysis view (no data copy)
        batchAna = buildBatchView(ana, batch_k_ok);

        % a. Blinks
        if RUN_BLINKS
            try; t0=tic; batchAna.computeBlinks(blinkOpts{:});
                fprintf('    [ok] blink (%.1fs)\n', toc(t0));
            catch e; fprintf('    [fail] blink: %s\n', e.message); end
        end

        % b. Signal stats
        if RUN_SIGNAL
            try; t0=tic; batchAna.computeSignalStats();
                fprintf('    [ok] signal (%.1fs)\n', toc(t0));
            catch e; fprintf('    [fail] signal: %s\n', e.message); end
        end

        % c. Quiet region
        if RUN_QUIET && ~isempty(batchAna.blink)
            try; t0=tic; batchAna.computeQuietRegion(quietOpts{:});
                fprintf('    [ok] quiet (%.1fs)\n', toc(t0));
            catch e; fprintf('    [fail] quiet: %s\n', e.message); end
        end

        % d. Spectral
        if RUN_SPECTRAL
            try; t0=tic; batchAna.computeSpectral(spectralOpts{:});
                fprintf('    [ok] spectral (%.1fs)\n', toc(t0));
            catch e; fprintf('    [fail] spectral: %s\n', e.message); end
        end

        % e. Temporal
        if RUN_TEMPORAL
            try; t0=tic; batchAna.computeTemporal(temporalOpts{:});
                fprintf('    [ok] temporal (%.1fs)\n', toc(t0));
            catch e; fprintf('    [fail] temporal: %s\n', e.message); end
        end

        % f. ICA
        if RUN_ICA
            try; t0=tic; batchAna.computeICA(icaOpts{:});
                fprintf('    [ok] ICA (%.1fs)\n', toc(t0));
            catch e; fprintf('    [fail] ICA: %s\n', e.message); end
        end

        % Harvest results into comboResults
        for bi = 1:numel(batch_k_ok)
            k = batch_k_ok(bi);
            comboResults(k) = harvestComboResult(batchAna, bi, ana.subject_results(k).parameters);
            doneFlags(k)    = true;
        end

        % Free signal memory for this batch
        if isSplit
            for k = batch_k_ok
                ana.unloadCombo(k);
            end
        end

        % Incremental save
        try
            save(saveFile, 'comboResults', 'doneFlags', '-v7.3');
            fprintf('    [saved] %d/%d combos done\n', sum(doneFlags), ana.nCombos);
        catch e
            fprintf('    [fail] save: %s\n', e.message);
        end

    end % batch loop

    % Build summary rows
    rows = buildSummaryRows(comboResults, subjectID, algoName, ana.nCombos, SUMMARY_METRICS);
    if ~isempty(rows), allSummaryRows{end+1} = rows; end %#ok<AGROW>

    logPrint(logFid, '  Total: %.1fs', toc(tSubject));
    fprintf('  Subject done (%.1fs)\n\n', toc(tSubject));

end % subject loop

%% ── WRITE SUMMARY TABLE ──────────────────────────────────────────────────

if ~isempty(allSummaryRows)
    try
        summaryTable = vertcat(allSummaryRows{:});
        writetable(summaryTable, fullfile(outputDir, 'summary_table.csv'));
        save(fullfile(outputDir, 'summary_table.mat'), 'summaryTable');
        fprintf('Summary table: %d rows x %d cols\n', height(summaryTable), width(summaryTable));
    catch e
        fprintf('[fail] Summary table: %s\n', e.message);
    end
end

fprintf('\n=== Batch complete in %.1fs ===\n  Output: %s\n\n', toc(tBatch), outputDir);
logPrint(logFid, '=== Batch complete: %.1fs ===', toc(tBatch));
fclose(logFid);


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function entries = discoverSubjects(dMeta, metaSuffix, dMono, monoSuffix)
    seen    = containers.Map('KeyType','char','ValueType','logical');
    entries = struct('subjectID',{},'algoName',{},'isSplit',{});
    for i = 1:numel(dMeta)
        base = dMeta(i).name(1:end-numel(metaSuffix));
        tok  = regexp(base, '^S(\d+)_(.+)$', 'tokens');
        if isempty(tok), continue; end
        sid = str2double(tok{1}{1}); algo = tok{1}{2};
        key = sprintf('%d|%s', sid, algo);
        if seen.isKey(key), continue; end
        seen(key) = true;
        entries(end+1).subjectID = sid; entries(end).algoName = algo; entries(end).isSplit = true; %#ok<AGROW>
    end
    for i = 1:numel(dMono)
        base = dMono(i).name(1:end-numel(monoSuffix));
        tok  = regexp(base, '^S(\d+)_(.+)$', 'tokens');
        if isempty(tok), continue; end
        sid = str2double(tok{1}{1}); algo = tok{1}{2};
        key = sprintf('%d|%s', sid, algo);
        if seen.isKey(key), continue; end
        seen(key) = true;
        entries(end+1).subjectID = sid; entries(end).algoName = algo; entries(end).isSplit = false; %#ok<AGROW>
    end
end

function cr = initComboResults(nCombos)
    cr = struct('parameters',cell(1,nCombos), ...
                'blink',cell(1,nCombos),'signalStats',cell(1,nCombos), ...
                'quietRegion',cell(1,nCombos),'spectral',cell(1,nCombos), ...
                'temporal',cell(1,nCombos),'ica',cell(1,nCombos));
end

function batchAna = buildBatchView(ana, batch_k)
    batchAna                  = ExperimentAnalysis(ana.accumulatedFolder, ...
                                    ana.algorithmName, ana.subjectID, ana.fs);
    batchAna.mode             = 'memory';
    batchAna.raw              = ana.raw;
    batchAna.subject_results  = ana.subject_results(batch_k);
    batchAna.nCombos          = numel(batch_k);
end

function cr = harvestComboResult(batchAna, bi, params)
    cr.parameters = params;
    for f = {'blink','signalStats','quietRegion','spectral','temporal','ica'}
        fn = f{1};
        obj = batchAna.(fn);
        if ~isempty(obj) && numel(obj.results) >= bi
            cr.(fn) = obj.results(bi);
        else
            cr.(fn) = [];
        end
    end
end

function rows = buildSummaryRows(comboResults, subjectID, algoName, nCombos, metrics)
    rows = table();
    if nCombos == 0, return; end
    rows.subjectID     = repmat(subjectID,  nCombos, 1);
    rows.algorithmName = repmat({algoName}, nCombos, 1);
    rows.comboIdx      = (1:nCombos)';
    if ~isempty(comboResults(1).parameters)
        pf = fieldnames(comboResults(1).parameters);
        for pi = 1:numel(pf)
            vals = NaN(nCombos,1);
            for c = 1:nCombos
                if ~isempty(comboResults(c).parameters)
                    v = comboResults(c).parameters.(pf{pi});
                    if isscalar(v) && isnumeric(v), vals(c) = v; end
                end
            end
            rows.(pf{pi}) = vals;
        end
    end
    for mi = 1:size(metrics,1)
        mName = metrics{mi,1}; fName = metrics{mi,2}; seg = metrics{mi,3};
        col   = [mName,'__',strrep(strrep(fName,'.','_'),' ','_')];
        if ~isempty(seg), col = [col,'__',seg]; end %#ok<AGROW>
        vals  = NaN(nCombos,1);
        for c = 1:nCombos
            try
                r = comboResults(c).(mName);
                if isempty(r), continue; end
                switch mName
                    case 'signalStats'
                        vals(c) = getNestedField(r, fName);
                    case {'blink','ica'}
                        vals(c) = getNestedField(r.summary, fName);
                    case {'spectral','temporal','quietRegion'}
                        if isempty(seg)
                            vals(c) = getNestedField(r.summary, fName);
                        else
                            vals(c) = getNestedField(r.summary.(seg), fName);
                        end
                end
            catch; end
        end
        rows.(col) = vals;
    end
end

function val = getNestedField(s, fieldPath)
    parts = strsplit(fieldPath, '.');
    val   = s;
    for i = 1:numel(parts), val = val.(parts{i}); end
end

function logPrint(fid, fmt, varargin)
    try; fprintf(fid,'[%s] %s\n',datestr(now,'HH:MM:SS'),sprintf(fmt,varargin{:})); catch; end %#ok<TNOW1,DATST>
end
