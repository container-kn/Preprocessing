%% emaASR — Parameter Sweep Template
% Runs a full factorial grid across multiple subjects.
% Results saved per-run to disk, then accumulated per subject, then cleaned up.
%
% Workflow:
%   configure -> sweep -> accumulate -> cleanup
%
% Analysis is done separately using ExperimentAnalysis on the accumulated folder.
%
% Key emaASR parameters to sweep:
%   cutoff          : SD threshold for artifact detection
%   beta1           : covariance matrix forgetting factor (lower = slower adaptation)
%   beta2           : threshold forgetting factor         (lower = slower adaptation)
%   cleanBuffLength : seconds of clean data buffer for baseline

clear; clc;

% ====================================================
% 1. Paths & Subjects
% ====================================================

projectRoot = 'J:\ASRProjct';
subjects    = [1 2];
fs          = 500;

% ====================================================
% 2. Base Configuration
% ====================================================

algorithm  = struct('name', 'ema-asr');
baseConfig = configExperiment(projectRoot, algorithm, 'fs', fs);

% ====================================================
% 3. Parameter Grid
% ====================================================
% Full factorial expansion — each field is a vector of values to sweep.
%
% Example: 2 cutoff x 2 beta1 = 4 combos x 2 subjects = 8 runs

% Note: blinks are typically 5-8 SD above baseline.
% cutoff=20 is too conservative to remove them — use 5-10 for blink reduction.
% Higher cutoffs (15-20) are appropriate if you only want to remove muscle/movement.
grid = struct( ...
    'cutoff', {[8:15]}, ...
    'beta1',  {[0.01:0.01:0.05]});

% beta2 typically matched to beta1 — to sweep independently:
% grid = struct( ...
%     'cutoff', {[10, 15, 20]}, ...
%     'beta1',  {[0.01, 0.03, 0.05]}, ...
%     'beta2',  {[0.01, 0.03, 0.05]});

% ====================================================
% 4. Run Sweep
% ====================================================

sweeper = parameterSweeper(baseConfig, algorithm, subjects);
sweeper.setSweep(grid);

% sweeper.sweeprun();                    % serial
sweeper.sweeprun('parallel', true); % parallel over combinations

% ====================================================
% 5. Accumulate Per Subject
% ====================================================
% Collects all per-run .mat files for each subject into one accumulated file.
%   results/accumulated/S1_ema-asr_accumulated.mat
%   results/accumulated/S2_ema-asr_accumulated.mat
%
% Each accumulated file contains:
%   raw             — raw signals stored ONCE (.calibration / .closed / .open)
%   subject_results — 1xN struct array, one entry per combo (cleaned only)

acc = accumulateResults( ...
    fullfile(projectRoot, 'results'), ...
    'ema-asr', subjects);

acc.reportEvery = 5;
acc.run();

% ====================================================
% 6. Cleanup
% ====================================================

acc.cleanup();

fprintf('\nSweep complete. Accumulated files in:\n  %s\n', acc.outputFolder);
