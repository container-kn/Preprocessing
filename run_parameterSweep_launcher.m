% run_parameterSweep_launcher.m
%
% Entry point called by run_asr_sweep.sh on PARAM Kamrupa.
%
% ── Parallelism design ───────────────────────────────────────────────────
% parameterSweeper.sweeprun() parallelises over combos only. With 10 subjects
% and 6 combos that gives 6 parallel tasks — wasting 42 of 48 cores.
%
% Instead, this launcher builds a flat list of ALL (subject, combo) pairs
% and runs a SINGLE parfor over all of them. With 10 subjects x 6 combos
% = 60 tasks, all 48 cores are immediately busy and the remaining 12 tasks
% queue behind them. No cores ever sit idle.
%
% ── Output files ─────────────────────────────────────────────────────────
% results/S<id>_<alg>_<params>.mat     one per (subject, combo) — from Experimenting
% results/S<id>_<alg>_raw.mat          raw signals, written once per subject
% results/accumulated/S<id>_<alg>.mat  merged per-subject file — from accumulateResults
%
% ── What to edit ─────────────────────────────────────────────────────────
%   Section 2 — ALGORITHM_NAME, SUBJECTS, signal settings
%   Section 3 — sweepGrid: parameters and value ranges to sweep
%   Section 1 — USE_PARALLEL: set false only when debugging locally

% =========================================================
% 0. Path setup
% =========================================================
PROJECT_ROOT = fullfile(getenv('HOME'), 'AdaptiveASR');

if ~isfolder(PROJECT_ROOT)
    error('launcher:badRoot', ...
        'PROJECT_ROOT not found: %s\nEdit this line to match your cluster path.', ...
        PROJECT_ROOT);
end

addpath(genpath(PROJECT_ROOT));

fprintf('\n=== ASR Sweep Launcher ===\n');
fprintf('  Project : %s\n', PROJECT_ROOT);
fprintf('  Started : %s\n\n', datestr(now));

% =========================================================
% 1. Parallelism — worker count from SLURM
% =========================================================
USE_PARALLEL = true;    % set false to run serially (debugging only)

nWorkers = str2double(getenv('SLURM_WORKERS'));
if isnan(nWorkers) || nWorkers < 1
    nWorkers = 1;       % running locally without SLURM
end

if USE_PARALLEL && nWorkers > 1
    fprintf('  Opening parpool: %d workers...\n', nWorkers);
    pool = parpool('local', nWorkers);
    fprintf('  Parpool ready.\n\n');
else
    USE_PARALLEL = false;
    fprintf('  Serial mode (nWorkers = %d).\n\n', nWorkers);
end

% =========================================================
% 2. Experiment configuration
% =========================================================
% ── Edit here ──
ALGORITHM_NAME = 'ema-asr';   % options: vanilla-asr | riemann-asr |
                               %          graph-asr  | hmo_asr | ema-asr

SUBJECTS = 1:10;               % subject IDs to process

algorithmStruct = struct('name', ALGORITHM_NAME);

baseConfig = configExperiment( ...
    PROJECT_ROOT,           ...
    algorithmStruct,        ...
    'fs',        500,       ...
    'bandpass',  [0.5 40],  ...
    'lineNoise', 50);

% =========================================================
% 3. Sweep grid
% =========================================================
% Each field = one parameter; value = vector of values to try.
% Full factorial: 3 cutoffs x 2 blocksizes = 6 combos.
%
% ema-asr params:   cutoff, blocksize, beta1, beta2, cleanBuffLength
% vanilla-asr:      cutoff, blocksize
% riemann-asr:      cutoff, blocksize, maxdims
% graph-asr:        cutoff, inflation, kNeighbors
% hmo_asr:          cutoff, stepsize

sweepGrid = struct( ...
    'cutoff',    [5, 8, 15], ...
    'blocksize', [10, 20]   ...
);

% =========================================================
% 4. Build full factorial task list: (subject, combo) pairs
% =========================================================
paramNames  = fieldnames(sweepGrid);
paramValues = struct2cell(sweepGrid);

% Unwrap any cell-of-cells
for i = 1:numel(paramValues)
    if iscell(paramValues{i})
        paramValues{i} = [paramValues{i}{:}];
    end
end

% Full factorial grid via ndgrid
grids = cell(1, numel(paramValues));
[grids{:}] = ndgrid(paramValues{:});
nCombos = numel(grids{1});

subjects = SUBJECTS(:)';
nSubj    = numel(subjects);
nTasks   = nSubj * nCombos;

% Flatten: each task is a (subjectIdx, comboIdx) pair
taskSubj  = zeros(1, nTasks);
taskCombo = zeros(1, nTasks);
t = 1;
for s = 1:nSubj
    for k = 1:nCombos
        taskSubj(t)  = subjects(s);
        taskCombo(t) = k;
        t = t + 1;
    end
end

fprintf('  Algorithm    : %s\n',  ALGORITHM_NAME);
fprintf('  Subjects     : %d\n',  nSubj);
fprintf('  Param combos : %d\n',  nCombos);
fprintf('  Total tasks  : %d\n',  nTasks);
fprintf('  Workers      : %d\n',  nWorkers);
fprintf('  Tasks/worker : %.1f\n\n', nTasks / max(nWorkers, 1));

% =========================================================
% 5. Run — parfor over ALL (subject, combo) pairs
% =========================================================
% Each worker independently:
%   1. Builds a configExperiment for its (subject, combo)
%   2. Runs Experimenting.run() which loads data, calibrates, processes
%   3. Saves results/S<id>_<alg>_<params>.mat
%   4. Writes results/S<id>_<alg>_raw.mat once (saveRawOnce guards duplication)
%
% Workers do NOT share state — each is fully independent.
% No broadcast variables — baseConfig_local is a plain struct copy.

baseConfig_local = baseConfig;   % copy for parfor broadcast (handle → struct copy)
algStruct_local  = algorithmStruct;
grids_local      = grids;
paramNames_local = paramNames;

if USE_PARALLEL

    parfor ti = 1:nTasks
        subjID    = taskSubj(ti);
        comboIdx  = taskCombo(ti);

        % Build parameter override for this combo
        override = struct();
        for i = 1:numel(paramNames_local)
            override.(paramNames_local{i}) = grids_local{i}(comboIdx);
        end

        fprintf('  [task %d/%d] S%d combo %d starting on worker %d\n', ...
            ti, nTasks, subjID, comboIdx, parallel.internal.pool.Communicator.worldRank);

        run_single_task(baseConfig_local, algStruct_local, subjID, override);
    end

else

    for ti = 1:nTasks
        subjID   = taskSubj(ti);
        comboIdx = taskCombo(ti);

        override = struct();
        for i = 1:numel(paramNames_local)
            override.(paramNames_local{i}) = grids_local{i}(comboIdx);
        end

        fprintf('  [%d/%d] S%d combo %d\n', ti, nTasks, subjID, comboIdx);
        run_single_task(baseConfig_local, algStruct_local, subjID, override);
    end

end

fprintf('\n=== All tasks complete: %s ===\n', datestr(now));

% =========================================================
% 6. Accumulate results into per-subject files
% =========================================================
fprintf('\n=== Accumulating results ===\n');

resultsFolder = fullfile(PROJECT_ROOT, 'results');
acc = accumulateResults(resultsFolder, ALGORITHM_NAME, subjects);
acc.run();

fprintf('=== Accumulation complete ===\n');

% =========================================================
% 7. Cleanup parpool
% =========================================================
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

fprintf('\n=== Launcher finished: %s ===\n', datestr(now));


% =========================================================
% Local helper — runs one (subject, combo) experiment
% =========================================================
function run_single_task(baseConfig, algStruct, subjID, override)
    % Wraps configExperiment + Experimenting.run() for one task.
    % Called from both serial and parfor loops.
    % All errors are caught and logged per-task so one failure doesn't
    % abort the whole sweep.

    try
        cfg = configExperiment( ...
            baseConfig.paths.project, ...
            algStruct,                ...
            'fs',              baseConfig.signal.fs,       ...
            'bandpass',        baseConfig.signal.bandpass, ...
            'lineNoise',       baseConfig.signal.lineNoise,...
            'subjectID',       subjID,                     ...
            'algorithmParams', override);

        cfg.sweepMode = true;   % raw saved once via saveRawOnce()

        exp = Experimenting(cfg);
        exp.run();              % load → preprocess → calibrate → process → save

    catch e
        % Print full error but don't rethrow — let other tasks continue
        fprintf('[ERROR] S%d override=%s\n  %s\n', ...
            subjID, ...
            strjoin(cellfun(@(f) sprintf('%s=%g', f, override.(f)), ...
                fieldnames(override), 'UniformOutput', false), ', '), ...
            getReport(e, 'extended'));
    end
end
