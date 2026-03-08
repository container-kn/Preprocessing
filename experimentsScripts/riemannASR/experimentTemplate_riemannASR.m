%% riemannASR — Single Subject Experiment Template
% Runs one subject with one parameter set end-to-end.
% Results saved to disk, blink counts shown as bar chart per channel,
% and zoomed blink windows plotted for visual inspection.
%
% Workflow:
%   configure -> run -> analyse in memory -> plots

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

algorithm        = struct('name', 'riemann-asr');
config           = configExperiment(projectRoot, algorithm, 'fs', fs);
config.subjectID = subjectID;

% Optional: override any default parameter
% config.algorithm.parameters.cutoff       = 15;
% config.algorithm.parameters.maxdims      = 0.66;  % fraction of dims to reconstruct
% config.algorithm.parameters.blocksize    = 10;
% config.algorithm.parameters.stepsize     = 32;    % processing step in samples
% config.algorithm.parameters.window_len   = 0.5;   % covariance window (sec)
% config.algorithm.parameters.window_overlap = 0.66;

fprintf('Config ready: Subject %d | Algorithm: %s\n', ...
    config.subjectID, config.algorithm.name);
fprintf('  cutoff=%d | maxdims=%.2f | blocksize=%d | stepsize=%d\n', ...
    config.algorithm.parameters.cutoff,    ...
    config.algorithm.parameters.maxdims,   ...
    config.algorithm.parameters.blocksize, ...
    config.algorithm.parameters.stepsize);

% ====================================================
% 3. Run Experiment
% ====================================================

experiment = Experimenting(config);
experiment.run();

% ====================================================
% 4. Extract Results
% ====================================================

[data, asrObj] = experiment.outputResults();

fprintf('Processed samples  : %d\n',   asrObj.nsamples);
fprintf('Modified fraction  : %.3f\n', mean(asrObj.modifiedMask));

% ====================================================
% 5. Blink Analysis
% ====================================================

ana = ExperimentAnalysis.fromOutputResults(data, asrObj, 'riemann-asr', fs);
ana.computeBlinks();

bSummary = ana.metrics(1).blink.summary;
fprintf('\n--- Blink Summary ---\n');
fprintf('  Raw blinks  (closed): %d  |  Clean: %d  |  Reduction: %.2f\n', ...
    bSummary.totalBlinks_raw_closed, bSummary.totalBlinks_clean_closed, bSummary.blinkReductionRatio_closed);
fprintf('  Raw blinks  (open)  : %d  |  Clean: %d  |  Reduction: %.2f\n', ...
    bSummary.totalBlinks_raw_open,   bSummary.totalBlinks_clean_open,   bSummary.blinkReductionRatio_open);

% ====================================================
% 6. Plot 1 — Blink counts per channel (bar chart)
% ====================================================

pc     = ana.metrics(1).blink.perChannel;
nChans = numel(pc.nBlinks_raw_closed);
chIdx  = 1:nChans;

figure('Name', sprintf('riemannASR S%d — Blink Counts per Channel', subjectID), ...
    'Position', [100 100 1000 420]);

subplot(1,2,1);
bar(chIdx, [pc.nBlinks_raw_closed(:), pc.nBlinks_clean_closed(:)]);
legend('Raw', 'Cleaned', 'Location', 'northeast');
title(sprintf('Closed Eyes  (total raw=%d  clean=%d)', ...
    bSummary.totalBlinks_raw_closed, bSummary.totalBlinks_clean_closed));
xlabel('Channel'); ylabel('Blink count');
xticks(chIdx); grid on;

subplot(1,2,2);
bar(chIdx, [pc.nBlinks_raw_open(:), pc.nBlinks_clean_open(:)]);
legend('Raw', 'Cleaned', 'Location', 'northeast');
title(sprintf('Open Eyes  (total raw=%d  clean=%d)', ...
    bSummary.totalBlinks_raw_open, bSummary.totalBlinks_clean_open));
xlabel('Channel'); ylabel('Blink count');
xticks(chIdx); grid on;

sgtitle(sprintf('riemannASR | Subject %d | cutoff=%d | maxdims=%.2f | blocksize=%d — Blink Counts', ...
    subjectID,                                 ...
    config.algorithm.parameters.cutoff,        ...
    config.algorithm.parameters.maxdims,       ...
    config.algorithm.parameters.blocksize));

% ====================================================
% 7. Plot 2 — Zoomed blink windows
% ====================================================

ch      = 1;
pad     = round(0.5 * fs);
nBlinks = 2;

blinkIdxClosed = pc.idx_raw_closed{ch};
blinkIdxOpen   = pc.idx_raw_open{ch};

figure('Name', sprintf('riemannASR S%d — Zoomed Blinks', subjectID), ...
    'Position', [100 560 1000 420]);

for seg = 1:2
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

    if isempty(bIdx)
        fprintf('No blinks detected in %s segment.\n', label);
        continue;
    end

    pick = round(linspace(1, numel(bIdx), min(nBlinks, numel(bIdx))));
    bIdx = bIdx(pick);

    for k = 1:numel(bIdx)
        iStart = max(1,               bIdx(k) - pad);
        iEnd   = min(size(raw,2),     bIdx(k) + pad);
        iEndC  = min(size(cleaned,2), bIdx(k) + pad);

        t_r = (iStart:iEnd)  / fs;
        t_c = (iStart:iEndC) / fs;

        subplot(2, nBlinks, (seg-1)*nBlinks + k);
        plot(t_r, raw(ch, iStart:iEnd),      'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
        plot(t_c, cleaned(ch, iStart:iEndC), 'b', 'LineWidth', 1.2);
        xline(bIdx(k)/fs, '--r', 'LineWidth', 0.8);
        title(sprintf('%s — blink %d', label, k));
        xlabel('Time (s)'); ylabel('Amplitude (\muV)');
        legend('Raw', 'Cleaned', 'Location', 'best');
        grid on;
    end
end

sgtitle(sprintf('riemannASR | Subject %d | cutoff=%d | maxdims=%.2f | blocksize=%d — Blink Windows (Ch %d)', ...
    subjectID,                                 ...
    config.algorithm.parameters.cutoff,        ...
    config.algorithm.parameters.maxdims,       ...
    config.algorithm.parameters.blocksize,     ...
    ch));
