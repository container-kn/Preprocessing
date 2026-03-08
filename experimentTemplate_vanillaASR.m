%% vanillaASR — Single Subject Experiment Template
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

projectRoot = 'J:\ArtifactSubspaceTracking';
% projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', 'ArtifactSubspaceTracking');

subjectID = 3;
fs        = 500;

% ====================================================
% 2. Configuration
% ====================================================

algorithm        = struct('name', 'vanilla-asr');
config           = configExperiment(projectRoot, algorithm, 'fs', fs);
config.subjectID = subjectID;

% Optional: override any default parameter
% config.algorithm.parameters.cutoff    = 15;
% config.algorithm.parameters.blocksize = 20;

fprintf('Config ready: Subject %d | Algorithm: %s | Cutoff: %d | Blocksize: %d\n', ...
    config.subjectID, ...
    config.algorithm.name, ...
    config.algorithm.parameters.cutoff, ...
    config.algorithm.parameters.blocksize);

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

ana = ExperimentAnalysis.fromOutputResults(data, asrObj, 'vanilla-asr', fs);
ana.computeBlinks();   % MAD-based bandpass detection (default)

b = ana.metrics(1).blink.summary;
fprintf('\n--- Blink Summary ---\n');
fprintf('  Raw blinks  (closed): %d  |  Clean: %d  |  Reduction: %.2f\n', ...
    b.totalBlinks_raw_closed, b.totalBlinks_clean_closed, b.blinkReductionRatio_closed);
fprintf('  Raw blinks  (open)  : %d  |  Clean: %d  |  Reduction: %.2f\n', ...
    b.totalBlinks_raw_open, b.totalBlinks_clean_open, b.blinkReductionRatio_open);

% ====================================================
% 6. Plot 1 — Blink counts per channel (bar chart)
% ====================================================

pc     = ana.metrics(1).blink.perChannel;
nChans = numel(pc.nBlinks_raw_closed);
chIdx  = 1:nChans;

figure('Name', sprintf('vanillaASR S%d — Blink Counts per Channel', subjectID), ...
    'Position', [100 100 1000 420]);

subplot(1,2,1);
bar(chIdx, [pc.nBlinks_raw_closed(:), pc.nBlinks_clean_closed(:)]);
legend('Raw', 'Cleaned', 'Location', 'northeast');
title(sprintf('Closed Eyes  (total raw=%d  clean=%d)', ...
    b.totalBlinks_raw_closed, b.totalBlinks_clean_closed));
xlabel('Channel'); ylabel('Blink count');
xticks(chIdx); grid on;

subplot(1,2,2);
bar(chIdx, [pc.nBlinks_raw_open(:), pc.nBlinks_clean_open(:)]);
legend('Raw', 'Cleaned', 'Location', 'northeast');
title(sprintf('Open Eyes  (total raw=%d  clean=%d)', ...
    b.totalBlinks_raw_open, b.totalBlinks_clean_open));
xlabel('Channel'); ylabel('Blink count');
xticks(chIdx); grid on;

sgtitle(sprintf('vanillaASR | Subject %d | Cutoff %d | Blocksize %d — Blink Counts', ...
    subjectID, ...
    config.algorithm.parameters.cutoff, ...
    config.algorithm.parameters.blocksize));

% ====================================================
% 7. Plot 2 — Zoomed blink windows
% ====================================================

ch      = 1;                  % channel to inspect
pad     = round(0.5 * fs);    % 0.5s either side of each blink peak
nBlinks = 2;                  % how many blink examples to show per segment

blinkIdxClosed = pc.idx_raw_closed{ch};
blinkIdxOpen   = pc.idx_raw_open{ch};

figure('Name', sprintf('vanillaASR S%d — Zoomed Blinks', subjectID), ...
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

    % Pick evenly-spaced examples across the segment
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

sgtitle(sprintf('vanillaASR | Subject %d | Cutoff %d | Blocksize %d — Blink Windows (Ch %d)', ...
    subjectID, ...
    config.algorithm.parameters.cutoff, ...
    config.algorithm.parameters.blocksize, ...
    ch));
