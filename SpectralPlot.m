classdef SpectralPlot < ASRPlotBase
    % SpectralPlot  PSD, band power, and scalp topography plots.
    %
    %   Inherits from ASRPlotBase.
    %   Reads from: ana.spectral  (SpectralAnalysis)
    %               ana.signalStats (SignalStatistics — for RRMSE/corr topo)
    %               ana.blink       (BlinkAnalysis    — for blink-reduction topo)
    %
    %   ── Plot methods ─────────────────────────────────────────────────────
    %
    %   plotPSD(c, seg)
    %       Overlaid raw vs clean mean PSD on a semilogy axis.
    %       Band regions (delta/theta/alpha/beta/gamma) shaded.
    %       Optionally plots per-channel PSDs as thin background lines.
    %
    %   plotPSDChannels(c, channels, seg)
    %       Small-multiples: one subplot per selected channel, each with
    %       raw vs clean PSD. Useful for identifying which channels change.
    %
    %   plotBandPower(c, seg)
    %       Grouped bar chart: raw vs clean absolute band power for each
    %       frequency band.  Change % labelled above each pair.
    %
    %   plotBandPowerChange(combos, seg)
    %       Line/bar: band power *change* ((clean-raw)/raw) per band,
    %       one line or group per combo. Shows which combos over-suppress.
    %
    %   plotTopoRRMSE(c, seg)          ** requires chanlocs + signalStats **
    %       Scalp topography of per-channel RRMSE (rms(clean-raw)/rms(raw)).
    %       Implements realScalpPlot.m approach for signal distortion metric.
    %
    %   plotTopoCorr(c, seg)           ** requires chanlocs + signalStats **
    %       Scalp topography of per-channel Pearson correlation raw vs clean.
    %
    %   plotTopoBlinkReduction(c, seg) ** requires chanlocs + blink **
    %       Scalp topography of per-channel blink count reduction fraction.
    %       Direct implementation of realScalpPlot.m.
    %
    %   plotTopoSweep(metric, seg)     ** requires chanlocs **
    %       Montage: one topoplot per parameter combo (subplot grid),
    %       replicating the nRows×nCols loop in realScalpPlot.m.
    %       metric = 'rrmse' | 'corr' | 'blinkReduction'
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       ana.computeBlinks();
    %       ana.computeSignalStats();
    %       sp = SpectralAnalysis(ana.raw, ana.subject_results, ana.fs, ...
    %                             ana.subjectID, ana.algorithmName);
    %       sp.compute();
    %       ana.spectral = sp;
    %
    %       tp = SpectralPlot(ana, EEG.chanlocs);
    %       tp.plotPSD(1, 'closed');
    %       tp.plotBandPower(1, 'closed');
    %       tp.plotTopoRRMSE(1, 'closed');
    %       tp.plotTopoBlinkReduction(1, 'closed');
    %       tp.plotTopoSweep('rrmse', 'closed');

    % ================================================================
    properties (Access = public)
        showPerChannelPSD = false  % show individual channel PSDs as thin lines
        bandNames = {'delta','theta','alpha','beta','gamma'}
        bandColors = [0.60 0.80 1.00;   % delta  — pale blue
                      0.80 1.00 0.80;   % theta  — pale green
                      1.00 1.00 0.65;   % alpha  — pale yellow
                      1.00 0.80 0.60;   % beta   — pale orange
                      1.00 0.70 0.70]   % gamma  — pale red
        topoStyle     = 'map'    % topoplot style: 'map' | 'contour' | 'both'
        topoColormap  = 'parula' % colormap for all topoplots
    end

    % ================================================================
    methods (Access = public)

        % ── Constructor ──────────────────────────────────────────────
        function obj = SpectralPlot(ana, chanlocs)
            if nargin < 2, chanlocs = []; end
            obj = obj@ASRPlotBase(ana, chanlocs);
        end

        % ════════════════════════════════════════════════════════════
        %  2. MEAN PSD  raw vs clean
        % ════════════════════════════════════════════════════════════
        function fig = plotPSD(obj, c, seg)
            % PLOTPSD  Mean PSD (raw vs clean) with shaded frequency bands.
            %
            %   c   : combo index
            %   seg : 'closed' | 'open' | 'all'

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end

            obj.hasSpectral(true);
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            psd = obj.ana.spectral.results(c).psd.(seg);
            f   = psd.f;
            Pr  = psd.raw_mean;
            Pc  = psd.clean_mean;

            % Per-channel PSDs for background (optional)
            if obj.showPerChannelPSD
                [Xraw, Xclean] = obj.getSegment(c, seg);
            end

            fig = obj.newFig(obj.makeTitle('PSD', seg, sprintf('combo %d', c)));
            ax  = axes(fig);
            hold(ax, 'on');

            % Band shading (behind everything)
            bandRanges = [1 4; 4 8; 8 13; 13 30; 30 100];
            fMax = min(f(end), 100);
            for bi = 1:5
                fl = bandRanges(bi,1);
                fh = min(bandRanges(bi,2), fMax);
                mask = f >= fl & f <= fh;
                if any(mask)
                    xl = [f(find(mask,1)), f(find(mask,1,'last'))];
                    patch(ax, [xl(1) xl(2) xl(2) xl(1)], ...
                          [1e-20 1e-20 1e20 1e20], ...
                          obj.bandColors(bi,:), ...
                          'EdgeColor','none','FaceAlpha',0.25, ...
                          'HandleVisibility','off');
                end
            end

            % Per-channel background lines
            if obj.showPerChannelPSD
                freqVec = 0 : obj.ana.spectral.freqResolution : obj.fs/2;
                nCh = size(Xraw, 1);
                for ch = 1:nCh
                    [PxxR, ~] = pwelch(Xraw(ch,:),   [], [], freqVec, obj.fs);
                    [PxxC, ~] = pwelch(Xclean(ch,:), [], [], freqVec, obj.fs);
                    semilogy(ax, freqVec, PxxR, '-', ...
                             'Color', [obj.colorRaw,   0.15], ...
                             'LineWidth', 0.5, 'HandleVisibility','off');
                    semilogy(ax, freqVec, PxxC, '-', ...
                             'Color', [obj.colorClean, 0.15], ...
                             'LineWidth', 0.5, 'HandleVisibility','off');
                end
            end

            % Mean PSDs
            semilogy(ax, f, Pr, '-', 'Color', obj.colorRaw,   ...
                     'LineWidth', obj.lineWidth+0.5, 'DisplayName', 'Raw');
            semilogy(ax, f, Pc, '--','Color', obj.colorClean,  ...
                     'LineWidth', obj.lineWidth+0.5, 'DisplayName', 'Clean');

            % Band name labels at top
            ax.XScale = 'log';
            bandMidFreq = [2 6 10.5 21 45];
            yl = [min([Pr;Pc])*0.5, max([Pr;Pc])*3];
            ylim(ax, yl);
            for bi = 1:5
                if bandMidFreq(bi) <= fMax
                    text(ax, bandMidFreq(bi), yl(2)*0.7, ...
                         obj.bandNames{bi}, ...
                         'HorizontalAlignment','center', ...
                         'FontSize', obj.fontSize-2, ...
                         'Color', [0.3 0.3 0.3]);
                end
            end

            xlim(ax, [1, fMax]);
            legend(ax, 'Location','southwest');
            obj.decorateAx(ax, 'Frequency (Hz)', 'Power (µV²/Hz)', ...
                obj.makeTitle('Mean PSD', seg, sprintf('combo %d', c)));
            obj.autoSave(fig, sprintf('psd_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  2. PSD PER CHANNEL — small multiples
        % ════════════════════════════════════════════════════════════
        function fig = plotPSDChannels(obj, c, channels, seg)
            % PLOTPSDCHANNELS  Per-channel PSDs as small-multiple subplots.
            %
            %   c        : combo index
            %   channels : channel indices
            %   seg      : 'closed' | 'open'

            if nargin < 3, channels = 1:min(9, size(obj.ana.raw.closed,1)); end
            if nargin < 4, seg = 'closed'; end

            obj.hasSpectral(true);
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            [Xraw, Xclean] = obj.getSegment(c, seg);
            C = size(Xraw, 1);
            obj.validateChannel(channels, C);

            nCh  = numel(channels);
            nC   = ceil(sqrt(nCh));
            nR   = ceil(nCh / nC);
            freqVec = 0 : obj.ana.spectral.freqResolution : obj.fs/2;
            fMax = min(freqVec(end), 100);

            fig = obj.newFig(obj.makeTitle( ...
                sprintf('PSD per channel (%d)', nCh), seg, ...
                sprintf('combo %d', c)));

            for ri = 1:nCh
                ch = channels(ri);
                ax = obj.subplotGrid(nR, nC, ri);
                hold(ax, 'on');

                [PxxR, f] = pwelch(double(Xraw(ch,:)),   [], [], freqVec, obj.fs);
                [PxxC, ~] = pwelch(double(Xclean(ch,:)), [], [], freqVec, obj.fs);

                semilogy(ax, f, PxxR, '-',  'Color', obj.colorRaw,   ...
                         'LineWidth', 1.0, 'DisplayName', 'Raw');
                semilogy(ax, f, PxxC, '--', 'Color', obj.colorClean,  ...
                         'LineWidth', 1.0, 'DisplayName', 'Clean');

                xlim(ax, [1, fMax]);
                title(ax, sprintf('Ch %d', ch), 'FontSize', obj.fontSize-1, ...
                      'Interpreter','none');
                ax.FontSize   = obj.fontSize - 2;
                ax.XTickLabel = {};
                ax.YTickLabel = {};
                grid(ax, 'on'); box(ax,'off');
                ax.XScale = 'log';

                if ri == 1
                    legend(ax, 'Location','southwest', 'FontSize', obj.fontSize-3);
                end
            end

            sgtitle(obj.makeTitle('PSD per channel', seg, ...
                sprintf('combo %d', c)), ...
                'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('psd_channels_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  2. BAND POWER — grouped bars
        % ════════════════════════════════════════════════════════════
        function fig = plotBandPower(obj, c, seg)
            % PLOTBANDPOWER  Grouped bar: raw vs clean per frequency band.
            %
            %   Shows absolute mean band power before and after ASR.
            %   Labels above each pair show % change.

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end

            obj.hasSpectral(true);
            obj.validateCombo(c);
            seg  = obj.resolveSegment(seg);

            pCh  = obj.ana.spectral.results(c).perChannel.(seg);
            bns  = obj.bandNames;
            nB   = numel(bns);

            rawPow   = zeros(nB, 1);
            cleanPow = zeros(nB, 1);
            for bi = 1:nB
                rawPow(bi)   = mean(pCh.raw.bandPower.(bns{bi}));
                cleanPow(bi) = mean(pCh.clean.bandPower.(bns{bi}));
            end

            fig = obj.newFig(obj.makeTitle( ...
                'Band power', seg, sprintf('combo %d', c)));
            ax  = axes(fig);

            bh = bar(ax, [rawPow, cleanPow]);
            bh(1).FaceColor = obj.colorRaw;
            bh(2).FaceColor = obj.colorClean;
            bh(1).FaceAlpha = 0.75;
            bh(2).FaceAlpha = 0.75;

            ax.XTick      = 1:nB;
            ax.XTickLabel = bns;

            % % change label above each pair
            xPos = bh(2).XEndPoints;
            yPos = max([rawPow, cleanPow], [], 2) * 1.06;
            for bi = 1:nB
                pct = 100 * (cleanPow(bi) - rawPow(bi)) / max(rawPow(bi), eps);
                col = obj.colorClean;
                if pct > 5,  col = [0.8 0 0]; end
                if pct < -5, col = [0 0.5 0]; end
                text(ax, xPos(bi), yPos(bi), sprintf('%+.0f%%', pct), ...
                     'HorizontalAlignment','center', ...
                     'FontSize', obj.fontSize-1, 'Color', col, ...
                     'FontWeight','bold');
            end

            legend(ax, {'Raw','Clean'}, 'Location','northeast');
            obj.decorateAx(ax, 'Frequency band', 'Mean power (µV²)', ...
                obj.makeTitle('Band power', seg, sprintf('combo %d', c)));
            obj.autoSave(fig, sprintf('bandpower_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  BAND POWER CHANGE — across combos
        % ════════════════════════════════════════════════════════════
        function fig = plotBandPowerChange(obj, combos, seg)
            % PLOTBANDPOWERCHANGE  (clean-raw)/raw per band for each combo.
            %
            %   combos : vector of combo indices (default = all)
            %   seg    : 'closed' | 'open'
            %
            %   Lines = bands, x-axis = combo index.

            if nargin < 2 || isempty(combos), combos = 1:obj.ana.nCombos; end
            if nargin < 3, seg = 'closed'; end

            obj.hasSpectral(true);
            seg  = obj.resolveSegment(seg);
            bns  = obj.bandNames;
            nB   = numel(bns);
            nCmb = numel(combos);

            changes = NaN(nB, nCmb);
            for ci = 1:nCmb
                c   = combos(ci);
                pCh = obj.ana.spectral.results(c).perChannel.(seg);
                for bi = 1:nB
                    rp = mean(pCh.raw.bandPower.(bns{bi}));
                    cp = mean(pCh.clean.bandPower.(bns{bi}));
                    changes(bi, ci) = (cp - rp) / max(rp, eps) * 100;
                end
            end

            fig = obj.newFig(obj.makeTitle('Band power change', seg));
            ax  = axes(fig);
            hold(ax, 'on');

            lineStyles = {'-o','-s','-^','-d','-v'};
            for bi = 1:nB
                plot(ax, 1:nCmb, changes(bi,:), lineStyles{bi}, ...
                     'DisplayName', bns{bi}, ...
                     'LineWidth', obj.lineWidth, ...
                     'MarkerSize', obj.markerSize);
            end

            yline(ax, 0, '--k', 'LineWidth', 0.8, 'HandleVisibility','off');
            ax.XTick = 1:nCmb;
            legend(ax, 'Location','best');
            obj.decorateAx(ax, 'Combo index', 'Power change (%)', ...
                obj.makeTitle('Band power change', seg));
            obj.autoSave(fig, sprintf('bandchange_%s', seg));
        end

        % ════════════════════════════════════════════════════════════
        %  3. TOPOPLOT — RRMSE  (requires signalStats + chanlocs)
        % ════════════════════════════════════════════════════════════
        function fig = plotTopoRRMSE(obj, c, seg)
            % PLOTTOPORRMSE  Scalp map of per-channel RRMSE.
            %
            %   RRMSE = rms(clean - raw) / rms(raw) per channel.
            %   High values = large distortion (bad).
            %   Source: ana.signalStats.results(c).segment.(seg).diff.rrmse_per_channel

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end

            obj.hasSignalStats(true);
            obj.checkTopoAvail();
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            vals = obj.ana.signalStats.results(c).segment.(seg).diff.rrmse_per_channel;
            vals = obj.ensureColumn(vals);

            fig = obj.topoFigure(vals, sprintf('RRMSE  |  %s  |  combo %d', seg, c), ...
                                 'RRMSE', [0, min(1, max(vals)*1.1)]);
            obj.autoSave(fig, sprintf('topo_rrmse_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  3. TOPOPLOT — CORRELATION
        % ════════════════════════════════════════════════════════════
        function fig = plotTopoCorr(obj, c, seg)
            % PLOTTOPOCORR  Scalp map of per-channel Pearson correlation.
            %
            %   corr = corr(raw, clean) per channel.
            %   Values close to 1 = well preserved signal (good).

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end

            obj.hasSignalStats(true);
            obj.checkTopoAvail();
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            vals = obj.ana.signalStats.results(c).segment.(seg).diff.corr_per_channel;
            vals = obj.ensureColumn(vals);

            fig = obj.topoFigure(vals, sprintf('Signal corr  |  %s  |  combo %d', seg, c), ...
                                 'Pearson r', [max(0, min(vals)-0.05), 1]);
            obj.autoSave(fig, sprintf('topo_corr_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  3. TOPOPLOT — BLINK REDUCTION  (realScalpPlot.m)
        % ════════════════════════════════════════════════════════════
        function fig = plotTopoBlinkReduction(obj, c, seg)
            % PLOTTOPOBLINKREDUCTION  Scalp map of per-channel blink reduction.
            %
            %   Per-channel reduction = (nRaw - nClean) / nRaw, clamped [0,1].
            %   Implements realScalpPlot.m: topoplot with parula colormap,
            %   clim [0 1], colorbar labelled 'Blink Reduction'.

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end

            obj.hasBlink(true);
            obj.checkTopoAvail();
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            perCh = obj.ana.blink.results(c).perChannel;
            nRaw  = double(perCh.(sprintf('nBlinks_raw_%s',   seg)));
            nCln  = double(perCh.(sprintf('nBlinks_clean_%s', seg)));

            % Per-channel fraction removed, clamped [0,1]
            reduction = zeros(numel(nRaw), 1);
            for ch = 1:numel(nRaw)
                if nRaw(ch) > 0
                    reduction(ch) = max(0, min(1, (nRaw(ch)-nCln(ch))/nRaw(ch)));
                end
            end

            fig = obj.topoFigure(reduction, ...
                sprintf('Blink reduction  |  %s  |  combo %d', seg, c), ...
                'Blink reduction', [0, 1]);
            obj.autoSave(fig, sprintf('topo_blinkreduction_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  3. TOPOPLOT SWEEP  — montage across combos
        %     Implements the nRows×nCols loop from realScalpPlot.m
        % ════════════════════════════════════════════════════════════
        function fig = plotTopoSweep(obj, metric, seg, combos)
            % PLOTTOPOSSWEEP  Montage of topoplots — one per combo.
            %
            %   metric : 'rrmse' | 'corr' | 'blinkReduction'
            %   seg    : 'closed' | 'open'
            %   combos : vector of combo indices (default = all)
            %
            %   Layout: ceil(sqrt(N)) × ceil(N/...) subplot grid.
            %   Shared colormap + colorbar on the right edge.
            %   Combo parameter value shown as each subplot title.
            %   Follows realScalpPlot.m design exactly.

            if nargin < 2, metric = 'rrmse';  end
            if nargin < 3, seg    = 'closed'; end
            if nargin < 4 || isempty(combos), combos = 1:obj.ana.nCombos; end

            obj.checkTopoAvail();
            seg  = obj.resolveSegment(seg);
            nCmb = numel(combos);

            % Collect per-channel vectors for each combo
            dataAll  = cell(nCmb, 1);
            titleStrs = cell(nCmb, 1);

            for ci = 1:nCmb
                c = combos(ci);
                dataAll{ci}   = obj.extractTopoMetric(metric, c, seg);
                titleStrs{ci} = obj.comboLabel(c);
            end

            % Shared color limits
            allVals = vertcat(dataAll{:});
            allVals = allVals(~isnan(allVals));
            switch metric
                case 'rrmse'
                    clim = [0, min(1, prctile(allVals, 95))];
                    cbLabel = 'RRMSE';
                case 'corr'
                    clim = [max(0, prctile(allVals, 5)), 1];
                    cbLabel = 'Pearson r';
                case 'blinkReduction'
                    clim = [0, 1];
                    cbLabel = 'Blink reduction';
                otherwise
                    clim = [min(allVals), max(allVals)];
                    cbLabel = metric;
            end

            nCols = ceil(sqrt(nCmb));
            nRows = ceil(nCmb / nCols);

            fig = obj.newFig(obj.makeTitle( ...
                sprintf('%s sweep', metric), seg));
            fig.Position(3) = nCols * 200;   % wider figure
            fig.Position(4) = nRows * 200;

            for ci = 1:nCmb
                ax = obj.subplotGrid(nRows, nCols, ci);
                vals = obj.ensureColumn(dataAll{ci});

                topoplot(vals, obj.chanlocs, 'style', obj.topoStyle, ...
                         'electrodes','off', 'headrad', 0.5);
                colormap(ax, obj.topoColormap);
                caxis(ax, clim);
                title(ax, titleStrs{ci}, ...
                      'FontSize', obj.fontSize-1, 'FontWeight','bold', ...
                      'Interpreter','none');
            end

            % Shared colorbar — positioned at right edge like realScalpPlot.m
            cb = colorbar('Position', [0.92 0.15 0.02 0.70]);
            cb.Label.String   = cbLabel;
            cb.Label.FontSize = obj.fontSize;
            colormap(obj.topoColormap);
            caxis(clim);

            sgtitle(obj.makeTitle(sprintf('%s sweep', metric), seg), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('topo_sweep_%s_%s', metric, seg));
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function vals = extractTopoMetric(obj, metric, c, seg)
            % EXTRACTTOPOMETRIC  Pull [C×1] per-channel values for given metric.
            switch lower(metric)
                case 'rrmse'
                    obj.hasSignalStats(true);
                    vals = obj.ana.signalStats.results(c).segment.(seg) ...
                               .diff.rrmse_per_channel;
                case 'corr'
                    obj.hasSignalStats(true);
                    vals = obj.ana.signalStats.results(c).segment.(seg) ...
                               .diff.corr_per_channel;
                case 'blinkreduction'
                    obj.hasBlink(true);
                    perCh = obj.ana.blink.results(c).perChannel;
                    nRaw  = double(perCh.(sprintf('nBlinks_raw_%s',   seg)));
                    nCln  = double(perCh.(sprintf('nBlinks_clean_%s', seg)));
                    vals  = zeros(numel(nRaw), 1);
                    for ch = 1:numel(nRaw)
                        if nRaw(ch) > 0
                            vals(ch) = max(0, min(1, (nRaw(ch)-nCln(ch))/nRaw(ch)));
                        end
                    end
                otherwise
                    error('SpectralPlot:unknownMetric', ...
                        'Unknown topo metric ''%s''. Use: rrmse, corr, blinkReduction.', ...
                        metric);
            end
            vals = obj.ensureColumn(vals);
        end

        function fig = topoFigure(obj, vals, titleStr, cbLabel, clim)
            % TOPOFIGURE  Single topoplot figure with colorbar.
            fig = obj.newFig(titleStr);
            ax  = axes(fig);

            topoplot(vals, obj.chanlocs, 'style', obj.topoStyle, ...
                     'electrodes','on', 'headrad', 0.5);
            colormap(ax, obj.topoColormap);

            if nargin >= 5 && ~isempty(clim)
                caxis(ax, clim);
            end

            obj.addTopoColorbar(ax, cbLabel, clim);
            title(ax, titleStr, 'FontSize', obj.fontSize, 'Interpreter','none');
        end

        function lbl = comboLabel(obj, c)
            % COMBOLABEL  Short parameter string for subplot title.
            params = obj.ana.subject_results(c).parameters;
            flds   = fieldnames(params);
            parts  = cell(numel(flds), 1);
            for i = 1:numel(flds)
                v = params.(flds{i});
                if isscalar(v) && isnumeric(v)
                    parts{i} = sprintf('%s=%.4g', flds{i}, v);
                end
            end
            parts = parts(~cellfun(@isempty, parts));
            if isempty(parts)
                lbl = sprintf('Combo %d', c);
            else
                lbl = strjoin(parts(1:min(2,end)), ' ');  % max 2 params
            end
        end

        function v = ensureColumn(~, v)
            v = v(:);
        end

    end % private methods

end
