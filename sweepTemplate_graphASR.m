%% graphASR — Incremental Parameter Sweep Template
% Runs a full factorial grid across multiple subjects and appends results
% into the existing split accumulated layout. Safe to re-run — already-
% accumulated combos are detected by parameter comparison and skipped.
%
% Workflow (every batch):
%   define grid -> run sweep -> patch (append new only) -> patchCleanup
%
% First run: patch() bootstraps the split layout from scratch.
% Subsequent runs: patch() appends only the new combos, skips duplicates.
%
% Analysis is done separately using ExperimentAnalysis on the accumulated folder.
%
% All graphASR parameters available to sweep:
%
%   GRAPH-SPECIFIC
%   inflation   : variance inflation threshold  (default 5)
%                 lower = more aggressive detection, higher = more conservative
%                 suggested range: [0.1, 0.25, 0.5, 0.74, 1, 1.5, 2, 2.5, 3, 4, 5, 8]
%   kNeighbors  : k-NN degree for channel graph construction  (default 4)
%                 controls which channels are considered neighbours
%                 suggested range: [1:6]
%   lambdaGain  : frequency-dependent detection gain multiplier  (default 5)
%                 scales spectral weighting in artifact scoring
%                 suggested range: [0.5, 1, 3, 5]
%   lambdaExp   : exponent shaping the frequency-gain curve  (default 0.5)
%                 low = flatter weighting, high = sharper roll-off
%                 suggested range: [0.25, 0.5, 1.0]
%   varKeep     : PCA variance retention for baseline covariance  (default 0.9999)
%                 lower drops noisy trailing components before calibration
%                 suggested range: [0.99, 0.999, 0.9999]
%   shrinkStyle : reconstruction mode  (default 'zero')
%                 'zero' = null out  |  'cap' = soft cap  |  'gamma' = regularised
%                 note: string field — use cell array syntax (see ADVANCED below)
%
%   SHARED BASE (rarely swept for graphASR, but available)
%   cutoff      : SD threshold for calibration distribution fit  (default 16)
%   blocksize   : covariance block size in samples  (default 20)

clear; clc;

% ====================================================
% 1. Paths & Subjects
% ====================================================

projectRoot = 'J:\ASRProjct';
subjects    = [6];
fs          = 500;

% ====================================================
% 2. Base Configuration
% ====================================================

algorithm  = struct('name', 'graph-asr');
baseConfig = configExperiment(projectRoot, algorithm, 'fs', fs);

% ====================================================
% 3. Parameter Grid  — edit this each batch
% ====================================================
% Each field is a vector of values for full-factorial expansion.
% Only combos NOT already in the accumulated layout will be run + patched.
%
% BATCH A — core 3-param sweep (12 × 6 × 4 = 288 combos)
% grid = struct( ...
%     'inflation',  {[0.1, 0.25, 0.5, 0.74, 1, 1.5, 2, 2.5, 3, 4, 5, 8]}, ...
%     'kNeighbors', {[1, 2, 3, 4, 5, 6]},                                   ...
%     'lambdaGain', {[0.5, 1, 3, 5]});

% BATCH B — add lambdaExp on top  (uncomment to use)
grid = struct( ...
    'inflation',  {[0.5, 1, 2, 3, 5, 8]}, ...    
    'lambdaGain', {[1, 3, 5]},             ...
    'lambdaExp',  {[0.25, 0.5, 1.0]});     % 6x3x3x3 = 162 combos

% BATCH C — varKeep sensitivity at your best point  (uncomment to use)
% grid = struct( ...
%     'inflation',  {[3]},                      ...
%     'kNeighbors', {[4]},                      ...
%     'lambdaGain', {[5]},                      ...
%     'varKeep',    {[0.99, 0.999, 0.9999]});

% ADVANCED — sweeping shrinkStyle (string field: must use cell array)
% Requires parameterSweeper to support cell-valued fields. Verify first.
% grid = struct( ...
%     'inflation',   {[1, 3, 5]},                  ...
%     'kNeighbors',  {[4]},                         ...
%     'lambdaGain',  {[3, 5]},                      ...
%     'shrinkStyle', {{'zero', 'cap', 'gamma'}});

% ====================================================
% 4. Run Sweep  — only new combos need running
% ====================================================

sweeper = parameterSweeper(baseConfig, algorithm, subjects);
sweeper.setSweep(grid);

% sweeper.sweeprun();                    % serial
sweeper.sweeprun('parallel', true);    % parallel over combinations

% ====================================================
% 5. Patch — append new combos into split accumulated layout
% ====================================================
% patch() reads the existing _meta.mat, skips any combo whose parameter set
% is already there, and appends only the genuinely new ones as new
% _combo_<k>.mat files. On first ever run it bootstraps the layout.
%
% Split layout written to:
%   results/accumulated/S<id>_graph-asr_raw.mat
%   results/accumulated/S<id>_graph-asr_meta.mat
%   results/accumulated/S<id>_graph-asr_combo_<k>.mat  (one per combo)

acc             = accumulateResults(fullfile(projectRoot, 'results'), 'graph-asr', subjects);
acc.splitFiles  = true;
acc.reportEvery = 10;

nNew = acc.patch();
fprintf('\n%d new combo(s) added.\n', nNew);

% Optional: print the full accumulated combo list to confirm
acc.listCombos();

% ====================================================
% 6. Cleanup — delete only the individual files just patched in
% ====================================================
% patchCleanup() only deletes files whose params are confirmed in _meta.mat.
% Files from a failed or partial run are left untouched.

acc.patchCleanup();

fprintf('\nBatch complete. Accumulated split files in:\n  %s\n', acc.outputFolder);
