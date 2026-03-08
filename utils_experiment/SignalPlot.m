classdef SignalPlot < ASRPlotBase
    % SignalPlot  Raw vs clean signal visualisation and blink event plots.
    %
    %   Inherits from ASRPlotBase (colour palette, figure helpers,
    %   data-access helpers, auto-save).
    %
    %   ── Plot methods ─────────────────────────────────────────────────────
    %
    %   plotWindow(c, channels, seg, timeWin)
    %       Multi-channel stacked view over any time window. Each channel
    %       gets its own subplot row: raw (blue) and clean (red), threshold
    %       dashed lines, and optional blink markers.  Implements and extends
    %       helperASR_plot_Signal to multiple channels.
    %
    %   plotBlinkZoom(c, ch, blinkNum, seg, windowSec)
    %       Single-channel zoomed view centred on one detected blink.
    %       raw, clean, threshold, and xline at blink centre.
    %       Implements helper_ASR_plot_BlinkEvent exactly.
    %
    %   plotBlinkOverlay(c, ch, seg, timeWin)
    %       Raw vs clean with blink peak markers scattered above the
    %       signal: raw blinks (black *), clean blinks (green *),
    %       persisted/matched blinks (orange o).
    %       Implements helperASR_plot_SignalANDblinks, extended with
    %       the matched-blink distinction from BlinkAnalysis.
    %
    %   plotBlinkCluster(c, ch, seg, windowSec)
    %       Finds the densest blink cluster in the recording and zooms
    %       to it.  Good for quickly finding the most interesting window.
    %       Annotates each blink with its status (removed / persisted).
    %
    %   plotAverageBlink(c, ch, seg, halfWinSec)
    %       Epoch all raw blinks ±halfWinSec, average across epochs.
    %       Overlay: mean raw (blue, ±1 std shading), mean clean (red).
    %       Zero-line at blink centre. Shows how ASR attenuates the average
    %       blink morphology.
    %
    %   plotMultiChannelStack(c, channels, seg, timeWin)
    %       Butterfly-style stacked plot: all selected channels in one axes
    %       with per-channel vertical offset.  Raw faint, clean solid.
    %       Good for overview — put adjacent to plotWindow for context.
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       ana = ExperimentAnalysis(folder, 'ema-asr', 3, 500);
    %       ana.load();
    %       ana.computeBlinks('method','mad');
    %
    %       sp = SignalPlot(ana);
    %       sp.saveDir = '/my/figs';           % optional auto-save
    %
    %       sp.plotWindow(1, [1 2 3], 'closed', [10 25]);
    %       sp.plotBlinkZoom(1, 1, 3, 'closed', 2);
    %       sp.plotBlinkOverlay(1, 1, 'closed', [0 60]);
    %       sp.plotBlinkCluster(1, 1, 'closed', 8);
    %       sp.plotAverageBlink(1, 1, 'closed', 0.5);
    %       sp.plotMultiChannelStack(1, 1:32, 'closed', [10 20]);
    %
    %       % With chanlocs passed at construction (enables topo in SpectralPlot)
    %       sp = SignalPlot(ana, EEG.chanlocs);

    % ================================================================
    properties (Access = public)
        showThreshold  = true   % overlay threshold line(s) on signal plots
        showBlinks     = true   % overlay blink markers when available
        channelOffset  = 150    % µV separation between channels in stack plot
        stackNormalize = false  % normalise each channel to unit std before stack
    end

    % ================================================================
    methods (Access = public)

        % ── Constructor ──────────────────────────────────────────────
        function obj = SignalPlot(ana, chanlocs)
            if nargin < 2, chanlocs = []; end
            obj = obj@ASRPlotBase(ana, chanlocs);
        end

        % ════════════════════════════════════════════════════════════
        %  1. MULTI-CHANNEL WINDOW PLOT
        %     helperASR_plot_Signal extended to N channels
        % ════════════════════════════════════════════════════════════
        function fig = plotWindow(obj, c, channels, seg, timeWin)
            % PLOTWINDOW  Stacked subplot per channel over a time window.
            %
            %   c        : combo index
            %   channels : channel indices, e.g. [1 2 3] or 1:32
            %   seg      : 'closed' | 'open' | 'calibration'
            %   timeWin  : [t_start t_end] seconds ([] = full segment)
            %
            %   Each row shows raw (blue) and clean (red).
            %   Threshold dashed line added if showThreshold=true and
            %   blink analysis has been run.
            %   Blink markers added if showBlinks=true.

            if nargin < 4, seg     = 'closed'; end
            if nargin < 5, timeWin = [];       end

            [Xraw, Xclean] = obj.getSegment(c, seg);
            [C, T]         = size(Xraw);
            obj.validateChannel(channels, C);

            t   = obj.timeAxis(T);
            if isempty(timeWin), timeWin = [t(1), t(end)]; end
            mask = obj.windowIdx(t, timeWin(1), timeWin(2));
            tw   = t(mask);

            thresh   = obj.getThreshold(c, seg);
            idxRaw   = obj.getBlinkIdx(c, seg, 'raw');
            idxClean = obj.getBlinkIdx(c, seg, 'clean');
            idxMatch = obj.getBlinkIdx(c, seg, 'matched_raw');

            nCh = numel(channels);
            fig = obj.newFig(obj.makeTitle( ...
                sprintf('%d ch', nCh), seg, ...
                sprintf('[%.1f–%.1f s]', timeWin(1), timeWin(2))));

            for ri = 1:nCh
                ch = channels(ri);
                ax = obj.subplotGrid(nCh, 1, ri);

                xr = Xraw(ch, mask);
                xc = Xclean(ch, mask);

                obj.drawSignalPair(ax, tw, xr, xc);

                % Threshold
                if obj.showThreshold && ~isempty(thresh)
                    obj.drawThreshold(ax, thresh(ch));
                end

                % Blink markers
                if obj.showBlinks
                    tBR = []; tBC = []; tBM = [];
                    if ~isempty(idxRaw)
                        ir = idxRaw{ch};
                        ir = ir(ir >= find(mask,1) & ir <= find(mask,1,'last'));
                        tBR = t(ir);
                    end
                    if ~isempty(idxClean)
                        ic = idxClean{ch};
                        ic = ic(ic >= find(mask,1) & ic <= find(mask,1,'last'));
                        tBC = t(ic);
                    end
                    if ~isempty(idxMatch)
                        im = idxMatch{ch};
                        im = im(im >= find(mask,1) & im <= find(mask,1,'last'));
                        tBM = t(im);
                    end
                    obj.drawBlinkMarkers(ax, tBR, tBC, tBM);
                end

                % Labels — only bottom axes gets x-label
                if ri == nCh
                    xLab = 'Time (s)';
                else
                    xLab = '';
                    ax.XTickLabel = {};
                end
                obj.decorateAx(ax, xLab, sprintf('Ch %d (µV)', ch), '');

                % Legend only on top axes
                if ri == 1
                    legend(ax, 'Location','northeast', 'FontSize', obj.fontSize-1);
                end
            end

            sgtitle(obj.makeTitle(seg, sprintf('combo %d', c)), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('window_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  1a. BLINK ZOOM  (helper_ASR_plot_BlinkEvent)
        % ════════════════════════════════════════════════════════════
        function fig = plotBlinkZoom(obj, c, ch, blinkNum, seg, windowSec)
            % PLOTBLINKZOOM  Zoomed single-channel view around one blink.
            %
            %   c         : combo index
            %   ch        : single channel index
            %   blinkNum  : which raw blink to centre on (1-indexed)
            %   seg       : 'closed' | 'open'
            %   windowSec : total width of window in seconds (default 2)
            %
            %   Shows: raw (blue), clean (red), threshold (--grey),
            %          xline at blink centre (--black).
            %   Equivalent to helper_ASR_plot_BlinkEvent.

            if nargin < 5, seg       = 'closed'; end
            if nargin < 6, windowSec = 2;        end

            obj.hasBlink(true);

            idxRaw = obj.getBlinkIdx(c, seg, 'raw');
            if isempty(idxRaw) || ch > numel(idxRaw) || isempty(idxRaw{ch})
                error('SignalPlot:noBlinksOnChannel', ...
                    'No raw blinks on channel %d segment ''%s''.', ch, seg);
            end

            blinkPeaks = idxRaw{ch};
            if blinkNum > numel(blinkPeaks)
                error('SignalPlot:blinkNumExceedsCount', ...
                    'blinkNum %d exceeds %d blinks on channel %d.', ...
                    blinkNum, numel(blinkPeaks), ch);
            end

            [Xraw, Xclean] = obj.getSegment(c, seg);
            T    = size(Xraw, 2);
            t    = obj.timeAxis(T);

            centerIdx  = blinkPeaks(blinkNum);
            halfSamp   = round((windowSec/2) * obj.fs);
            idxRange   = max(1, centerIdx-halfSamp) : min(T, centerIdx+halfSamp);
            tw         = t(idxRange);

            thresh = obj.getThreshold(c, seg);

            fig = obj.newFig(obj.makeTitle( ...
                sprintf('Blink zoom  ch%d  blink #%d', ch, blinkNum), seg));

            ax = axes(fig);
            obj.drawSignalPair(ax, tw, Xraw(ch,idxRange), Xclean(ch,idxRange));

            if obj.showThreshold && ~isempty(thresh)
                obj.drawThreshold(ax, thresh(ch));
            end

            % Centre line
            xline(ax, t(centerIdx), '--k', 'LineWidth', 1.2, ...
                  'DisplayName', 'Blink centre');

            % Mark blink number annotation
            yl = ylim(ax);
            text(ax, t(centerIdx), yl(2)*0.95, sprintf('Blink #%d', blinkNum), ...
                 'HorizontalAlignment','center', 'FontSize', obj.fontSize-1, ...
                 'Color', obj.colorBlink);

            legend(ax, 'Location','northeast');
            obj.decorateAx(ax, 'Time (s)', 'Amplitude (µV)', ...
                obj.makeTitle(sprintf('Blink zoom  ch%d  blink #%d', ch, blinkNum), ...
                              seg, sprintf('combo %d', c)));
            obj.autoSave(fig, sprintf('blinkzoom_c%d_ch%d_b%d_%s', c, ch, blinkNum, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  1a. BLINK OVERLAY  (helperASR_plot_SignalANDblinks)
        % ════════════════════════════════════════════════════════════
        function fig = plotBlinkOverlay(obj, c, ch, seg, timeWin)
            % PLOTBLINKOVERLAY  Raw vs clean with blink scatter markers.
            %
            %   c       : combo index
            %   ch      : single channel index
            %   seg     : 'closed' | 'open'
            %   timeWin : [t_start t_end] in seconds ([] = full segment)
            %
            %   Raw blinks    — black  *  (peaks in original signal)
            %   Clean blinks  — green  *  (peaks still present after ASR)
            %   Persisted     — orange o  (matched: blink that survived)
            %   Equivalent to helperASR_plot_SignalANDblinks, extended
            %   to distinguish persisted vs newly-created clean blinks.

            if nargin < 4, seg     = 'closed'; end
            if nargin < 5, timeWin = [];       end

            obj.hasBlink(true);

            [Xraw, Xclean] = obj.getSegment(c, seg);
            [C, T]         = size(Xraw);
            obj.validateChannel(ch, C);

            t    = obj.timeAxis(T);
            if isempty(timeWin), timeWin = [t(1), t(end)]; end
            mask = obj.windowIdx(t, timeWin(1), timeWin(2));
            tw   = t(mask);

            % Collect blink times in window
            tRaw = obj.blinkTimesInWindow(c, seg, 'raw',          ch, timeWin, t);
            tCln = obj.blinkTimesInWindow(c, seg, 'clean',        ch, timeWin, t);
            tMat = obj.blinkTimesInWindow(c, seg, 'matched_raw',  ch, timeWin, t);

            thresh = obj.getThreshold(c, seg);

            fig = obj.newFig(obj.makeTitle( ...
                sprintf('Blink overlay  ch%d', ch), seg, ...
                sprintf('[%.1f–%.1f s]', timeWin(1), timeWin(2))));

            ax = axes(fig);
            obj.drawSignalPair(ax, tw, Xraw(ch,mask), Xclean(ch,mask));

            if obj.showThreshold && ~isempty(thresh)
                obj.drawThreshold(ax, thresh(ch));
            end

            obj.drawBlinkMarkers(ax, tRaw, tCln, tMat);

            legend(ax, 'Location','northeast');
            obj.decorateAx(ax, 'Time (s)', 'Amplitude (µV)', ...
                obj.makeTitle(sprintf('Blink overlay  ch%d', ch), ...
                              seg, sprintf('combo %d', c)));
            obj.autoSave(fig, sprintf('blinkoverlay_c%d_ch%d_%s', c, ch, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  1a. BLINK CLUSTER  — find densest window automatically
        % ════════════════════════════════════════════════════════════
        function fig = plotBlinkCluster(obj, c, ch, seg, windowSec)
            % PLOTBLINKCLUSTER  Auto-find the densest blink cluster and plot it.
            %
            %   Slides a window of windowSec seconds across the raw blink
            %   times, finds the window containing the most blinks, then
            %   calls plotBlinkOverlay on that window.
            %   Each blink annotated: 'R' removed, 'P' persisted, 'N' new.
            %
            %   c         : combo index
            %   ch        : single channel
            %   seg       : 'closed' | 'open'
            %   windowSec : cluster window width (default 8 s)

            if nargin < 4, seg       = 'closed'; end
            if nargin < 5, windowSec = 8;        end

            obj.hasBlink(true);

            idxRaw = obj.getBlinkIdx(c, seg, 'raw');
            if isempty(idxRaw) || ch > numel(idxRaw) || isempty(idxRaw{ch})
                error('SignalPlot:noBlinksOnChannel', ...
                    'No raw blinks on channel %d segment ''%s''.', ch, seg);
            end

            [Xraw, ~] = obj.getSegment(c, seg);
            T     = size(Xraw, 2);
            t_all = obj.timeAxis(T);

            tRawAll = idxRaw{ch} / obj.fs;

            % Sliding window count
            halfW   = windowSec / 2;
            bestCount = 0;
            bestCentre = median(tRawAll);
            for i = 1:numel(tRawAll)
                tc  = tRawAll(i);
                cnt = sum(tRawAll >= tc-halfW & tRawAll <= tc+halfW);
                if cnt > bestCount
                    bestCount  = cnt;
                    bestCentre = tc;
                end
            end

            t0 = max(t_all(1),  bestCentre - halfW);
            t1 = min(t_all(end), bestCentre + halfW);

            % Re-use plotBlinkOverlay
            fig = obj.plotBlinkOverlay(c, ch, seg, [t0, t1]);

            % Add annotation labels per blink
            ax     = fig.CurrentAxes;
            idxCln = obj.getBlinkIdx(c, seg, 'clean');
            idxMat = obj.getBlinkIdx(c, seg, 'matched_raw');

            tC   = []; if ~isempty(idxCln), tC = idxCln{ch}/obj.fs; end
            tM   = []; if ~isempty(idxMat), tM = idxMat{ch}/obj.fs; end
            tInW = tRawAll(tRawAll >= t0 & tRawAll <= t1);

            yl  = ylim(ax);
            ylo = yl(1) + 0.05*(yl(2)-yl(1));
            for i = 1:numel(tInW)
                tb = tInW(i);
                if any(abs(tM - tb) < 0.1)
                    lbl = 'P';  col = obj.colorMatched;
                elseif ~isempty(tC) && ~any(abs(tC - tb) < 0.3)
                    lbl = 'R';  col = obj.colorBlink;
                else
                    lbl = 'P';  col = obj.colorMatched;
                end
                text(ax, tb, ylo, lbl, 'Color', col, ...
                     'FontSize', obj.fontSize, 'FontWeight','bold', ...
                     'HorizontalAlignment','center');
            end

            % Update title
            ax.Title.String = obj.makeTitle( ...
                sprintf('Blink cluster  ch%d  (%d blinks)', ch, bestCount), ...
                seg, sprintf('combo %d', c));
            obj.autoSave(fig, sprintf('blinkcluster_c%d_ch%d_%s', c, ch, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  5. AVERAGE BLINK WAVEFORM
        % ════════════════════════════════════════════════════════════
        function fig = plotAverageBlink(obj, c, ch, seg, halfWinSec)
            % PLOTAVERAGEBLINK  Average blink morphology across all raw blinks.
            %
            %   Epochs every raw blink ±halfWinSec and averages across events.
            %   Shows:
            %     Raw   mean (blue, solid)   ± std (blue shading)
            %     Clean mean (red,  solid)   ± std (red  shading)
            %     Zero-line at blink centre
            %     Dashed threshold at ±threshold level
            %
            %   This reveals whether ASR consistently zeros the blink peak
            %   (good) or leaves a residual (bad) or distorts the baseline.
            %
            %   c          : combo index
            %   ch         : single channel
            %   seg        : 'closed' | 'open'
            %   halfWinSec : epoch half-width in seconds (default 0.5)

            if nargin < 4, seg        = 'closed'; end
            if nargin < 5, halfWinSec = 0.5;      end

            obj.hasBlink(true);

            idxRaw = obj.getBlinkIdx(c, seg, 'raw');
            if isempty(idxRaw) || ch > numel(idxRaw) || isempty(idxRaw{ch})
                error('SignalPlot:noBlinks', ...
                    'No raw blinks on channel %d, segment ''%s''.', ch, seg);
            end

            [Xraw, Xclean] = obj.getSegment(c, seg);
            T       = size(Xraw, 2);
            halfN   = round(halfWinSec * obj.fs);
            tEpoch  = (-halfN:halfN) / obj.fs;
            nSamp   = numel(tEpoch);

            peaks = idxRaw{ch};

            % Build epoch matrices [nBlinks × nSamples]
            epochsRaw   = NaN(numel(peaks), nSamp);
            epochsClean = NaN(numel(peaks), nSamp);

            for ei = 1:numel(peaks)
                i0 = peaks(ei) - halfN;
                i1 = peaks(ei) + halfN;
                if i0 < 1 || i1 > T, continue; end
                epochsRaw(ei,:)   = Xraw(ch,   i0:i1);
                epochsClean(ei,:) = Xclean(ch, i0:i1);
            end

            % Remove epochs that were clipped (all-NaN rows)
            valid       = ~any(isnan(epochsRaw), 2);
            epochsRaw   = epochsRaw(valid,:);
            epochsClean = epochsClean(valid,:);
            nValid      = sum(valid);

            muR  = mean(epochsRaw,   1);
            sigR = std(epochsRaw,    0, 1);
            muC  = mean(epochsClean, 1);
            sigC = std(epochsClean,  0, 1);

            thresh = obj.getThreshold(c, seg);

            fig = obj.newFig(obj.makeTitle( ...
                sprintf('Average blink  ch%d  (n=%d)', ch, nValid), seg));
            ax = axes(fig);
            hold(ax, 'on');

            % Shading first (under lines)
            obj.shadeBand(ax, tEpoch, muR, sigR);
            c2 = obj.colorClean * 0.6 + [0.4 0.4 0.4];  % lighter red
            fill(ax, [tEpoch, fliplr(tEpoch)], ...
                     [muC+sigC, fliplr(muC-sigC)], c2, ...
                     'EdgeColor','none','FaceAlpha',0.30, ...
                     'HandleVisibility','off');

            % Mean lines
            plot(ax, tEpoch, muR, '-', 'Color', obj.colorRaw,   ...
                 'LineWidth', obj.lineWidth+0.5, 'DisplayName', 'Raw mean');
            plot(ax, tEpoch, muC, '-', 'Color', obj.colorClean, ...
                 'LineWidth', obj.lineWidth+0.5, 'DisplayName', 'Clean mean');

            % Centre line
            xline(ax, 0, '--k', 'LineWidth', 1, 'DisplayName', 'Blink centre');

            % Threshold
            if obj.showThreshold && ~isempty(thresh)
                yline(ax,  thresh(ch), '--', 'Color', obj.colorThresh, ...
                      'LineWidth', 0.9, 'DisplayName', 'Threshold');
                yline(ax, -thresh(ch), '--', 'Color', obj.colorThresh, ...
                      'LineWidth', 0.9, 'HandleVisibility','off');
            end

            legend(ax, 'Location','northeast');
            obj.decorateAx(ax, 'Time from blink centre (s)', 'Amplitude (µV)', ...
                obj.makeTitle(sprintf('Average blink  ch%d  n=%d', ch, nValid), ...
                              seg, sprintf('combo %d', c)));
            obj.autoSave(fig, sprintf('avgblink_c%d_ch%d_%s', c, ch, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  MULTI-CHANNEL STACK (butterfly overlay with offsets)
        % ════════════════════════════════════════════════════════════
        function fig = plotMultiChannelStack(obj, c, channels, seg, timeWin)
            % PLOTMULTICHANNELSTACK  All channels in a single axes, offset.
            %
            %   Each channel is shifted vertically by channelOffset µV so
            %   traces don't overlap. Raw is faint, clean is solid.
            %   Good for a quick whole-montage overview.
            %
            %   c         : combo index
            %   channels  : channel indices (e.g. 1:32)
            %   seg       : 'closed' | 'open'
            %   timeWin   : [t_start t_end] seconds ([] = full)

            if nargin < 4, seg     = 'closed'; end
            if nargin < 5, timeWin = [];       end

            [Xraw, Xclean] = obj.getSegment(c, seg);
            [C, T]         = size(Xraw);
            obj.validateChannel(channels, C);

            t    = obj.timeAxis(T);
            if isempty(timeWin), timeWin = [t(1), t(end)]; end
            mask = obj.windowIdx(t, timeWin(1), timeWin(2));
            tw   = t(mask);

            nCh    = numel(channels);
            offset = obj.channelOffset;

            fig = obj.newFig(obj.makeTitle( ...
                'Channel stack', seg, ...
                sprintf('[%.1f–%.1f s]', timeWin(1), timeWin(2))));
            ax = axes(fig);
            hold(ax, 'on');

            yTicks = zeros(nCh, 1);
            for ri = 1:nCh
                ch  = channels(ri);
                off = (ri-1) * offset;
                xr  = Xraw(ch,   mask);
                xc  = Xclean(ch, mask);

                if obj.stackNormalize
                    s  = max(std(xr), eps);
                    xr = xr / s * (offset*0.45);
                    xc = xc / s * (offset*0.45);
                end

                plot(ax, tw, xr + off, '-', ...
                     'Color', [obj.colorRaw, 0.35], ...
                     'LineWidth', 0.8, 'HandleVisibility','off');
                plot(ax, tw, xc + off, '-', ...
                     'Color', obj.colorClean, ...
                     'LineWidth', 1.0, 'HandleVisibility','off');
                yTicks(ri) = off;
            end

            % Dummy lines for legend
            plot(ax, NaN, NaN, '-', 'Color', obj.colorRaw,   ...
                 'LineWidth', 1.2, 'DisplayName', 'Raw');
            plot(ax, NaN, NaN, '-', 'Color', obj.colorClean, ...
                 'LineWidth', 1.2, 'DisplayName', 'Clean');

            ax.YTick      = yTicks;
            ax.YTickLabel = arrayfun(@(ch) sprintf('Ch%d', ch), channels, ...
                                     'UniformOutput', false);
            legend(ax, 'Location','northeast');
            obj.decorateAx(ax, 'Time (s)', '', ...
                obj.makeTitle('Channel stack', seg, sprintf('combo %d', c)));
            obj.autoSave(fig, sprintf('chstack_c%d_%s', c, seg));
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function tBlinks = blinkTimesInWindow(obj, c, seg, role, ch, timeWin, t)
            % BLINKTIMESINWINDOW  Blink times (s) on channel ch within timeWin.
            tBlinks = [];
            idx = obj.getBlinkIdx(c, seg, role);
            if isempty(idx) || ch > numel(idx) || isempty(idx{ch})
                return;
            end
            peaks   = idx{ch};
            % Clip to valid sample range
            peaks   = peaks(peaks >= 1 & peaks <= numel(t));
            times   = t(peaks);
            tBlinks = times(times >= timeWin(1) & times <= timeWin(2));
        end

    end % private methods

end
