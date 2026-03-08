%% graphASR — Single Subject Experiment Template
% Runs one subject with one parameter set end-to-end.
% Results saved to disk and optionally analysed in-memory immediately.
%
% Workflow:
%   configure -> run -> (optional) analyse in memory
%
% To sweep parameters across subjects instead, see:
%   graphASR_sweep_template.m

clear; clc;

% ====================================================
% 1. Paths & Identity
% ====================================================

projectRoot = 'J:\ASRProjct';
subjectID   = 3;
fs          = 500;

% ====================================================
% 2. Configuration
% ====================================================

algorithm        = struct('name', 'graph-asr');
config           = configExperiment(projectRoot, algorithm, 'fs', fs);
config.subjectID = subjectID;

% Optional: override any default parameter
% config.algorithm.parameters.cutoff     = 15;
% config.algorithm.parameters.kNeighbors = 6;
% config.algorithm.parameters.lambdaGain = 5;
% config.algorithm.parameters.inflation  = 5;
% config.algorithm.parameters.shrinkStyle = 'cap';   % 'zero' | 'cap' | 'gamma'

fprintf('Config ready: Subject %d | Algorithm: %s\n', ...
    config.subjectID, config.algorithm.name);
fprintf('  cutoff=%d | kNeighbors=%d | lambdaGain=%d | inflation=%d\n', ...
    config.algorithm.parameters.cutoff,     ...
    config.algorithm.parameters.kNeighbors, ...
    config.algorithm.parameters.lambdaGain, ...
    config.algorithm.parameters.inflation);

% ====================================================
% 3. Run Experiment
% ====================================================

experiment = Experimenting(config);
experiment.run();

% ====================================================
% 4. Extract Results (in-memory)
% ====================================================

[data, asrObj] = experiment.outputResults();

fprintf('Processed samples  : %d\n',   asrObj.nsamples);
fprintf('Modified fraction  : %.3f\n', mean(asrObj.modifiedMask));
fprintf('Graph nodes (chans): %d\n',   asrObj.nchans);
fprintf('Retained PCA rank  : %d\n',   asrObj.r);

% ====================================================
% 5. Quick In-Memory Analysis (optional)
% ====================================================

ana = ExperimentAnalysis.fromOutputResults(data, asrObj, 'graph-asr', fs);
ana.computeBlinks();

b = ana.metrics(1).blink.summary;
fprintf('\n--- Blink Summary ---\n');
fprintf('  Raw blinks  (closed): %d  |  Clean: %d  |  Reduction: %.2f\n', ...
    b.totalBlinks_raw_closed, b.totalBlinks_clean_closed, b.blinkReductionRatio_closed);
fprintf('  Raw blinks  (open)  : %d  |  Clean: %d  |  Reduction: %.2f\n', ...
    b.totalBlinks_raw_open,   b.totalBlinks_clean_open,   b.blinkReductionRatio_open);

% ====================================================
% 6. Graph-specific probe inspection (optional)
% ====================================================

if ~isempty(asrObj.probeRaw)
    pr  = asrObj.probeRaw;
    nW  = pr.step;     % number of windows recorded

    fprintf('\n--- Graph Probe (Raw) ---\n');
    fprintf('  Windows recorded      : %d\n',   nW);
    fprintf('  Mean spectral entropy : %.3f\n', mean(pr.spectralEntropy(1:nW), 'omitnan'));
    fprintf('  Mean Riemann drift    : %.3f\n', mean(pr.riemannDrift(1:nW),    'omitnan'));
    fprintf('  Mean subspace angle   : %.2f deg\n', ...
        mean(pr.subspaceAngle(1:nW), 'omitnan'));
    fprintf('  Mean spectral disagree: %.3f\n', mean(pr.spectralDisagreement(1:nW), 'omitnan'));
end

% ====================================================
% 7. Quick Plots — blink-centred windows only
% ====================================================

ch      = 1;       % channel to inspect
pad     = round(0.5 * fs);   % 0.5s either side of each blink peak
nBlinks = 2;       % how many blinks to show per segment

% Pull blink indices from the analysis probe (raw signal, closed eyes)
blinkIdxClosed = ana.metrics(1).blink.perChannel.idx_raw_closed{ch};
blinkIdxOpen   = ana.metrics(1).blink.perChannel.idx_raw_open{ch};

figure('Name', sprintf('graphASR — S%d (blink windows)', subjectID));

for seg = 1:2    % 1 = closed, 2 = open
    if seg == 1
        raw     = data.closedSignal;
        cleaned = data.cleanedClosed;
        bIdx    = blinkIdxClosed;
        label   = 'Closed Eyes';
    else
        raw     = data.openSignal;
        cleaned = data.cleanedOpen;
        bIdx    = blinkIdxOpen;
        label   = 'Open Eyes';
    end

    % Pick up to nBlinks evenly spaced blink events
    if isempty(bIdx)
        fprintf('No blinks detected in %s segment', label);
        continue;
    end
    pick = round(linspace(1, numel(bIdx), min(nBlinks, numel(bIdx))));
    bIdx = bIdx(pick);

    for b = 1:numel(bIdx)
        iStart = max(1,                  bIdx(b) - pad);
        iEnd   = min(size(raw, 2),       bIdx(b) + pad);
        iEndC  = min(size(cleaned, 2),   bIdx(b) + pad);

        t_r = (iStart:iEnd)   / fs;
        t_c = (iStart:iEndC)  / fs;

        subplot(2, nBlinks, (seg-1)*nBlinks + b);
        plot(t_r, raw(ch, iStart:iEnd),         'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
        plot(t_c, cleaned(ch, iStart:iEndC),    'b',  'LineWidth', 1.2);
        xline(bIdx(b)/fs, '--r', 'LineWidth', 0.8);
        title(sprintf('%s — blink %d', label, b));
        xlabel('Time (s)'); ylabel('Amplitude');
        legend('Raw', 'Cleaned', 'Location', 'best');
        grid on;
    end
end

sgtitle(sprintf('graphASR | Subject %d | cutoff=%d | kN=%d | lGain=%d | infl=%d', ...
    subjectID,                              ...
    config.algorithm.parameters.cutoff,     ...
    config.algorithm.parameters.kNeighbors, ...
    config.algorithm.parameters.lambdaGain, ...
    config.algorithm.parameters.inflation));
