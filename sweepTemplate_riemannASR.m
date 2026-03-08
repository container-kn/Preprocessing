%% riemannASR — Parameter Sweep Template
% Runs a full factorial grid across multiple subjects.
% Results saved per-run to disk, then accumulated per subject, then cleaned up.
%
% Workflow:
%   configure -> sweep -> accumulate -> cleanup
%
% Analysis is done separately using ExperimentAnalysis on the accumulated folder.
%
% Key riemannASR parameters to sweep:
%   cutoff  : SD threshold — same role as in vanillaASR
%   maxdims : fraction of subspace dimensions to reconstruct (0 < maxdims <= 1)
%             lower = more aggressive dimensionality reduction

clear; clc;

% ====================================================
% 1. Paths & Subjects
% ====================================================

projectRoot = 'J:\ASRProjct';
subjects    = [1];
fs          = 500;

% ====================================================
% 2. Base Configuration
% ====================================================

algorithm  = struct('name', 'riemann-asr');
baseConfig = configExperiment(projectRoot, algorithm, 'fs', fs);

% ====================================================
% 3. Parameter Grid
% ====================================================
% Full factorial expansion — each field is a vector of values to sweep.
%
% Example: 3 cutoff x 2 maxdims = 6 combos x 2 subjects = 12 runs

grid = struct( ...
    'cutoff',  {[10, 20]}, ...
    'maxdims', {[0.5, 0.66]});

% To also sweep blocksize:
% grid = struct( ...
%     'cutoff',    {[10, 15, 20]}, ...
%     'maxdims',   {[0.5, 0.66]},  ...
%     'blocksize', {[10, 20]});    % 3x2x2 = 12 combos

% ====================================================
% 4. Run Sweep
% ====================================================

sweeper = parameterSweeper(baseConfig, algorithm, subjects);
sweeper.setSweep(grid);

sweeper.sweeprun();                    % serial
% sweeper.sweeprun('parallel', true); % parallel over combinations

% ====================================================
% 5. Accumulate Per Subject
% ====================================================
% Collects all per-run .mat files for each subject into one accumulated file.
%   results/accumulated/S1_riemann-asr_accumulated.mat
%   results/accumulated/S2_riemann-asr_accumulated.mat
%
% Each accumulated file contains:
%   raw             — raw signals stored ONCE (.calibration / .closed / .open)
%   subject_results — 1xN struct array, one entry per combo (cleaned only)

acc = accumulateResults( ...
    fullfile(projectRoot, 'results'), ...
    'riemann-asr', subjects);

acc.reportEvery = 5;
acc.run();

% ====================================================
% 6. Cleanup
% ====================================================

acc.cleanup();

fprintf('\nSweep complete. Accumulated files in:\n  %s\n', acc.outputFolder);
