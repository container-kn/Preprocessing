%% vanillaASR — Parameter Sweep Template
% Runs a full factorial grid of cutoff x blocksize across multiple subjects.
% Results saved per-run to disk, then accumulated per subject, then cleaned up.
%
% Workflow:
%   configure -> sweep -> accumulate -> cleanup
%
% Analysis is done separately using ExperimentAnalysis on the accumulated folder.

clear; clc;

% ====================================================
% 1. Paths & Subjects
% ====================================================

projectRoot = 'J:\ASRProjct';
subjects    = [6:10];
fs          = 500;

% ====================================================
% 2. Base Configuration
% ====================================================

algorithm  = struct('name', 'vanilla-asr');
baseConfig = configExperiment(projectRoot, algorithm, 'fs', fs);

% ====================================================
% 3. Parameter Grid
% ====================================================
% Full factorial expansion — each field is a vector of values to sweep.
%   cutoff    : SD threshold (lower = more aggressive cleaning)
%   blocksize : covariance block size in samples

grid = struct( ...
    'cutoff',    {[0:2:6,7:1:15,18:2:25,30:10:100]});

% 2 cutoff x 1 blocksize = 2 combos x 2 subjects = 4 total runs

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
%   results/accumulated/S1_vanilla-asr_accumulated.mat
%   results/accumulated/S2_vanilla-asr_accumulated.mat
%
% Each accumulated file contains:
%   raw             — raw signals stored ONCE (.calibration / .closed / .open)
%   subject_results — 1xN struct array, one entry per combo (cleaned only)

acc = accumulateResults( ...
    fullfile(projectRoot, 'results'), ...
    'vanilla-asr', subjects);
acc.splitFiles = true;


acc.reportEvery = 5;
acc.run();
acc.resplit()
% ====================================================
% 6. Cleanup
% ====================================================
% Deletes individual per-run .mat files after confirming each accumulated
% file was written successfully.

acc.cleanup();

fprintf('\nSweep complete. Accumulated files in:\n  %s\n', acc.outputFolder);
