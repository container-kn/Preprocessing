%% graphASR — Parameter Sweep Template
% Runs a full factorial grid across multiple subjects.
% Results saved per-run to disk, then accumulated per subject, then cleaned up.
%
% Workflow:
%   configure -> sweep -> accumulate -> cleanup
%
% Analysis is done separately using ExperimentAnalysis on the accumulated folder.
%
% Key graphASR parameters to sweep:
%   inflation   : variance inflation threshold — controls detection sensitivity
%                 (lower = more aggressive, higher = more conservative)
%   kNeighbors  : k-NN degree for channel graph construction
%   lambdaGain  : frequency-dependent detection gain
%   shrinkStyle : 'zero' | 'cap' | 'gamma'  (usually fixed, not swept)

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
% 3. Parameter Grid
% ====================================================
% Full factorial expansion — each field is a vector of values to sweep.
%
% Example below: 3 inflation x 2 kNeighbors = 6 combos x 2 subjects = 12 runs

grid = struct( ...
    'inflation',  {[0.1, 0.25, 0.5, 0.74, 1, 1.5,2, 2.5, 3, 4, 5, 8]}, ...
    'kNeighbors', {[1:6]}, ...
    'lambdaGain', {[0.5, 1, 3, 5]});

% To also sweep lambdaGain:
% grid = struct( ...
%     'inflation',  {[3, 5, 8]}, ...
%     'kNeighbors', {[4, 6]},    ...
%     'lambdaGain', {[3, 5]});   % 3x2x2 = 12 combos

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
%   results/accumulated/S1_graph-asr_accumulated.mat
%   results/accumulated/S2_graph-asr_accumulated.mat
%
% Each accumulated file contains:
%   raw             — raw signals stored ONCE (.calibration / .closed / .open)
%   subject_results — 1xN struct array, one entry per combo (cleaned only)

acc = accumulateResults( ...
    fullfile(projectRoot, 'results'), ...
    'graph-asr', subjects);

acc.reportEvery = 5;
acc.run();

% ====================================================
% 6. Cleanup
% ====================================================
% Deletes individual per-run .mat files after confirming each accumulated
% file was written successfully.

acc.cleanup();

fprintf('\nSweep complete. Accumulated files in:\n  %s\n', acc.outputFolder);
