classdef SummaryPlot < ASRPlotBase
    % SummaryPlot  Parameter-sweep and cross-combo summary visualisation.
    %
    %   Inherits from ASRPlotBase.
    %   Reads from: ana.blink / ana.signalStats / ana.spectral / ana.temporal
    %
    %   Designed for answering: "How does metric X change as I sweep
    %   parameter Y across all combos?"
    %
    %   ── Plot methods ─────────────────────────────────────────────────────
    %
    %   plotBlinkSweep(seg)
    %       4-panel: blink reduction %, persistence rate, removal rate,
    %       artefact rate — all vs combo index. Best single figure for
    %       blink summary across a parameter sweep.
    %
    %   plotSignalQualitySweep(seg)
    %       3-panel: RRMSE, Pearson corr, SNR(dB) vs combo index.
    %
    %   plotMetricSweep(analyser, fieldName, seg)
    %       Generic single-metric line plot vs combo index.
    %       Works for any getSummary()-compatible field.
    %
    %   plotParameterGrid(xField, metrics, seg)
    %       Scatter: x-axis = one scalar parameter (e.g. 'cutoff'),
    %       y-axis = metric value.  One subplot per metric.
    %       Useful when combos are a 1-D sweep over one parameter.
    %
    %   plotComboHeatmap(metrics, seg)
    %       Heatmap: combos × metrics, z = normalised metric value.
    %       Colour = parula, rows = combos, columns = metrics.
    %       Good for comparing many metrics at once.
    %
    %   plotTemporalSweep(fieldName, seg)
    %       Temporal metric (getSummary from TemporalAnalysis) vs combo.
    %       Adds reference line at 0 and labels each point.
    %
    %   plotBandPowerSweep(band, seg)
    %       Power change in one frequency band vs combo index.
    %       Shows raw power and clean power side by side.
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       sp = SummaryPlot(ana);
    %       sp.saveDir = '/my/figs';
    %
    %       sp.plotBlinkSweep('closed');
    %       sp.plotSignalQualitySweep('closed');
    %       sp.plotMetricSweep('blink', 'blinkReductionPct_closed', '');
    %       sp.plotParameterGrid('cutoff', ...
    %           {{'blink','blinkReductionPct_closed',''}, ...
    %            {'signalStats','rrmse','closed','diff'}}, 'closed');
    %       sp.plotComboHeatmap(...);
    %       sp.plotTemporalSweep('energyCompaction.dimReduction_mean', 'closed');

    % ================================================================
    properties (Access = public)
        showComboLabels = true   % annotate each point with combo label
        maxLabelChars   = 14     % truncate long combo labels
        comboParamField = ''     % if set, use this field for x-axis tick labels
                                 % e.g. 'cutoff'; '' = use combo index
    end

    % ================================================================
    methods (Access = public)

        % ── Constructor ──────────────────────────────────────────────
        function obj = SummaryPlot(ana, chanlocs)
            if nargin < 2, chanlocs = []; end
            obj = obj@ASRPlotBase(ana, chanlocs);
        end

        % ════════════════════════════════════════════════════════════
        %  BLINK SWEEP — 4-panel
        % ════════════════════════════════════════════════════════════
        function fig = plotBlinkSweep(obj, seg)
            % PLOTBLINKSWEEP  4-panel blink metrics across combos.
            %
            %   Panels: reduction %, persistence rate,
            %           removal rate, artefact rate.
            %   seg : 'closed' | 'open' (default: 'closed')

            if nargin < 2, seg = 'closed'; end
            obj.hasBlink(true);
            seg = obj.resolveSegment(seg);

            metrics = { ...
                sprintf('blinkReductionPct_%s',  seg), 'Reduction (%)';
                sprintf('persistenceRate_%s',     seg), 'Persistence rate';
                sprintf('removalRate_%s',         seg), 'Removal rate';
                sprintf('artefactRate_%s',        seg), 'Artefact rate'};

            refLines = {NaN; 0; 1; 0};   % ideal target values

            fig = obj.newFig(obj.makeTitle('Blink sweep', seg));
            x   = obj.comboXAxis();

            for pi = 1:4
                vals = obj.ana.blink.getSummary(metrics{pi,1});
                ax   = obj.subplotGrid(2, 2, pi);
                obj.sweepLine(ax, x, vals, metrics{pi,2}, refLines{pi});
                obj.setComboXTicks(ax, x);
                obj.decorateAx(ax, 'Combo', metrics{pi,2}, ...
                    sprintf('%s  |  %s', metrics{pi,2}, seg));
            end

            sgtitle(obj.makeTitle('Blink metrics sweep', seg), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('blinksweep_%s', seg));
        end

        % ════════════════════════════════════════════════════════════
        %  SIGNAL QUALITY SWEEP — 3-panel
        % ════════════════════════════════════════════════════════════
        function fig = plotSignalQualitySweep(obj, seg)
            % PLOTSIGNALQUALITYSWEEP  RRMSE, correlation, SNR across combos.
            %
            %   seg : 'closed' | 'open'

            if nargin < 2, seg = 'closed'; end
            obj.hasSignalStats(true);
            seg = obj.resolveSegment(seg);

            metrics = { ...
                'rrmse',        'diff',  'RRMSE',   0;
                'corr',         'diff',  'Corr (r)', 1;
                'snrDB',        'diff',  'SNR (dB)', NaN};

            fig = obj.newFig(obj.makeTitle('Signal quality sweep', seg));
            x   = obj.comboXAxis();

            for pi = 1:3
                vals = obj.ana.signalStats.getSummary( ...
                    metrics{pi,1}, seg, metrics{pi,2});
                ax = obj.subplotGrid(1, 3, pi);
                obj.sweepLine(ax, x, vals, metrics{pi,3}, metrics{pi,4});
                obj.setComboXTicks(ax, x);
                obj.decorateAx(ax, 'Combo', metrics{pi,3}, ...
                    sprintf('%s  |  %s', metrics{pi,3}, seg));
            end

            sgtitle(obj.makeTitle('Signal quality sweep', seg), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('signalquality_%s', seg));
        end

        % ════════════════════════════════════════════════════════════
        %  GENERIC SINGLE-METRIC SWEEP
        % ════════════════════════════════════════════════════════════
        function fig = plotMetricSweep(obj, analyser, fieldName, seg, signalType)
            % PLOTMETRICSWEEP  Any getSummary-compatible field vs combo.
            %
            %   analyser   : 'blink' | 'signalStats' | 'spectral' | 'temporal'
            %   fieldName  : field name / dot-path passed to getSummary()
            %   seg        : segment ('' or [] for analysers that don't need it)
            %   signalType : (signalStats only) 'raw'|'clean'|'diff'

            if nargin < 4, seg        = 'closed'; end
            if nargin < 5, signalType = 'diff';   end

            vals = obj.pullMetric(analyser, fieldName, seg, signalType);
            x    = obj.comboXAxis();

            label = strrep(fieldName, '_', ' ');
            fig   = obj.newFig(obj.makeTitle(label, seg));
            ax    = axes(fig);

            obj.sweepLine(ax, x, vals, label, NaN);
            obj.setComboXTicks(ax, x);
            obj.decorateAx(ax, 'Combo', label, ...
                obj.makeTitle(label, seg));
            obj.autoSave(fig, sprintf('sweep_%s_%s_%s', analyser, ...
                strrep(fieldName,'.','_'), seg));
        end

        % ════════════════════════════════════════════════════════════
        %  PARAMETER GRID — scatter vs one scalar parameter
        % ════════════════════════════════════════════════════════════
        function fig = plotParameterGrid(obj, xField, metricSpecs, seg)
            % PLOTPARAMETERGRID  Metrics plotted against one swept parameter.
            %
            %   xField      : parameter name in subject_results.parameters,
            %                 e.g. 'cutoff' or 'blocksize'
            %   metricSpecs : cell array of {analyser, fieldName, seg, signalType}
            %                 e.g. {{'blink','blinkReductionPct_closed','',''}, ...
            %                       {'signalStats','rrmse','closed','diff'}}
            %   seg         : default segment if not specified per metric

            if nargin < 4, seg = 'closed'; end

            % Extract x-axis values
            xVals = NaN(obj.ana.nCombos, 1);
            for c = 1:obj.ana.nCombos
                p = obj.ana.subject_results(c).parameters;
                if isfield(p, xField)
                    v = p.(xField);
                    if isscalar(v) && isnumeric(v)
                        xVals(c) = v;
                    end
                end
            end

            if all(isnan(xVals))
                error('SummaryPlot:noParam', ...
                    'Parameter ''%s'' not found in subject_results.parameters.', ...
                    xField);
            end

            nM  = numel(metricSpecs);
            nC  = ceil(sqrt(nM));
            nR  = ceil(nM / nC);
            fig = obj.newFig(obj.makeTitle( ...
                sprintf('vs %s', xField), seg));

            for mi = 1:nM
                spec    = metricSpecs{mi};
                mAna    = spec{1};
                mField  = spec{2};
                mSeg    = spec{3};  if isempty(mSeg),  mSeg  = seg;    end
                mSigTyp = '';
                if numel(spec) >= 4, mSigTyp = spec{4}; end

                vals = obj.pullMetric(mAna, mField, mSeg, mSigTyp);
                ax   = obj.subplotGrid(nR, nC, mi);
                hold(ax, 'on');

                scatter(ax, xVals, vals, obj.markerSize*8, ...
                        obj.colorRaw, 'filled', 'HandleVisibility','off');

                % Fit a line if enough points
                validMask = ~isnan(xVals) & ~isnan(vals(:)');
                if sum(validMask) >= 3
                    p = polyfit(xVals(validMask), vals(validMask), 1);
                    xFit = linspace(min(xVals(validMask)), ...
                                    max(xVals(validMask)), 50);
                    plot(ax, xFit, polyval(p,xFit), '-', ...
                         'Color', obj.colorClean, 'LineWidth', 1.2, ...
                         'HandleVisibility','off');
                end

                label = strrep(mField,'_',' ');
                obj.decorateAx(ax, xField, label, ...
                    sprintf('%s  |  %s', label, mSeg));
            end

            sgtitle(obj.makeTitle(sprintf('Metrics vs %s', xField), seg), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('paramgrid_%s_%s', xField, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  COMBO HEATMAP
        % ════════════════════════════════════════════════════════════
        function fig = plotComboHeatmap(obj, analyserMetricList, seg)
            % PLOTCOMBOHEATMAP  combos × metrics heatmap, z = normed value.
            %
            %   analyserMetricList : cell of {analyser, fieldName, [seg], [sigType]}
            %     e.g. { {'blink','blinkReductionPct_closed'},
            %             {'signalStats','rrmse','closed','diff'},
            %             {'temporal','signal.rrmse_mean','closed'} }
            %   seg : default segment

            if nargin < 3, seg = 'closed'; end

            nM = numel(analyserMetricList);
            nC = obj.ana.nCombos;

            data   = NaN(nC, nM);
            labels = cell(1, nM);

            for mi = 1:nM
                spec   = analyserMetricList{mi};
                mAna   = spec{1};
                mField = spec{2};
                mSeg   = seg;
                mSig   = 'diff';
                if numel(spec) >= 3, mSeg = spec{3}; end
                if numel(spec) >= 4, mSig = spec{4}; end

                try
                    v = obj.pullMetric(mAna, mField, mSeg, mSig);
                    data(:, mi) = v(:);
                catch
                    % leave NaN
                end

                labels{mi} = obj.shortLabel(mField);
            end

            % Normalise columns to [0,1] for display
            dataNorm = data;
            for mi = 1:nM
                col = data(:, mi);
                lo  = min(col); hi = max(col);
                if hi > lo
                    dataNorm(:, mi) = (col - lo) / (hi - lo);
                end
            end

            fig = obj.newFig(obj.makeTitle('Combo heatmap', seg));
            ax  = axes(fig);

            imagesc(ax, dataNorm);
            colormap(ax, 'parula');
            cb = colorbar(ax);
            cb.Label.String = 'Normalised value';

            ax.XTick      = 1:nM;
            ax.XTickLabel = labels;
            ax.XTickLabelRotation = 35;
            ax.YTick      = 1:nC;
            ax.YTickLabel = arrayfun(@(c) obj.comboShortLabel(c), ...
                                     1:nC, 'UniformOutput', false);
            ax.FontSize   = obj.fontSize - 1;

            % Annotate cells with actual values
            for ci = 1:nC
                for mi = 1:nM
                    v = data(ci, mi);
                    if ~isnan(v)
                        text(ax, mi, ci, sprintf('%.2g', v), ...
                             'HorizontalAlignment','center', ...
                             'VerticalAlignment','middle', ...
                             'FontSize', obj.fontSize-3, 'Color','k');
                    end
                end
            end

            obj.decorateAx(ax, 'Metric', 'Combo', ...
                obj.makeTitle('Combo heatmap', seg));
            obj.autoSave(fig, sprintf('heatmap_%s', seg));
        end

        % ════════════════════════════════════════════════════════════
        %  TEMPORAL METRIC SWEEP
        % ════════════════════════════════════════════════════════════
        function fig = plotTemporalSweep(obj, fieldName, seg)
            % PLOTTEMPORALSWEEP  TemporalAnalysis summary field vs combo.
            %
            %   fieldName : dot-delimited getSummary path, e.g.
            %       'signal.rrmse_mean'
            %       'energyCompaction.dimReduction_mean'
            %       'calibrationDrift.energyCaptureGain_mean'
            %   seg : 'closed' | 'open'

            if nargin < 3, seg = 'closed'; end
            obj.hasTemporal(true);
            seg = obj.resolveSegment(seg);

            vals  = obj.ana.temporal.getSummary(fieldName, seg);
            x     = obj.comboXAxis();
            label = strrep(fieldName, '.', '  ›  ');
            label = strrep(label, '_', ' ');

            fig = obj.newFig(obj.makeTitle('Temporal sweep', fieldName, seg));
            ax  = axes(fig);
            hold(ax, 'on');

            obj.sweepLine(ax, x, vals, label, 0);

            % Annotate improvement direction
            if contains(fieldName,'Gain') || contains(fieldName,'dimReduction')
                ylabel_str = sprintf('%s\n(↑ better)', label);
            elseif contains(fieldName,'rrmse') || contains(fieldName,'angle')
                ylabel_str = sprintf('%s\n(↓ better)', label);
            else
                ylabel_str = label;
            end

            obj.setComboXTicks(ax, x);
            obj.decorateAx(ax, 'Combo', ylabel_str, ...
                obj.makeTitle('Temporal sweep', fieldName, seg));
            obj.autoSave(fig, sprintf('temporal_%s_%s', ...
                strrep(fieldName,'.','_'), seg));
        end

        % ════════════════════════════════════════════════════════════
        %  BAND POWER SWEEP — one frequency band
        % ════════════════════════════════════════════════════════════
        function fig = plotBandPowerSweep(obj, band, seg)
            % PLOTBANDPOWERSWEEP  Raw vs clean mean power in one band.
            %
            %   band : 'delta' | 'theta' | 'alpha' | 'beta' | 'gamma'
            %   seg  : 'closed' | 'open'

            if nargin < 2, band = 'alpha'; end
            if nargin < 3, seg  = 'closed'; end

            obj.hasSpectral(true);
            seg = obj.resolveSegment(seg);
            x   = obj.comboXAxis();

            rawPow   = NaN(1, obj.ana.nCombos);
            cleanPow = NaN(1, obj.ana.nCombos);
            for c = 1:obj.ana.nCombos
                pCh = obj.ana.spectral.results(c).perChannel.(seg);
                rawPow(c)   = mean(pCh.raw.bandPower.(band));
                cleanPow(c) = mean(pCh.clean.bandPower.(band));
            end

            fig = obj.newFig(obj.makeTitle( ...
                sprintf('%s band power sweep', band), seg));
            ax  = axes(fig);
            hold(ax, 'on');

            plot(ax, x, rawPow,   '-o', 'Color', obj.colorRaw,   ...
                 'LineWidth', obj.lineWidth, 'MarkerSize', obj.markerSize, ...
                 'DisplayName', 'Raw');
            plot(ax, x, cleanPow, '-s', 'Color', obj.colorClean,  ...
                 'LineWidth', obj.lineWidth, 'MarkerSize', obj.markerSize, ...
                 'DisplayName', 'Clean');

            legend(ax, 'Location','best');
            obj.setComboXTicks(ax, x);
            obj.decorateAx(ax, 'Combo', sprintf('%s power (µV²)', band), ...
                obj.makeTitle(sprintf('%s band power', band), seg));
            obj.autoSave(fig, sprintf('bandpowersweep_%s_%s', band, seg));
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        % ── Core sweep line plot with optional reference line ─────

        function sweepLine(obj, ax, x, vals, label, refVal)
            % SWEEPLINE  Line + markers + optional horizontal reference.

            hold(ax, 'on');

            % Fill NaN gaps with visible gap in line
            plot(ax, x, vals, '-o', ...
                 'Color',     obj.colorRaw, ...
                 'LineWidth', obj.lineWidth, ...
                 'MarkerFaceColor', obj.colorRaw, ...
                 'MarkerSize', obj.markerSize, ...
                 'HandleVisibility','off');

            % Reference line
            if nargin >= 6 && ~isnan(refVal)
                yline(ax, refVal, '--', 'Color', [0.5 0.5 0.5], ...
                      'LineWidth', 0.9, ...
                      'Label', sprintf('ideal=%.3g', refVal), ...
                      'LabelHorizontalAlignment','left', ...
                      'HandleVisibility','off');
            end

            % Per-point labels
            if obj.showComboLabels
                for ci = 1:numel(x)
                    lbl = obj.comboShortLabel(ci);
                    text(ax, x(ci), vals(ci), lbl, ...
                         'VerticalAlignment','bottom', ...
                         'HorizontalAlignment','center', ...
                         'FontSize', obj.fontSize-3, ...
                         'Color', [0.35 0.35 0.35], ...
                         'Interpreter','none');
                end
            end
        end

        % ── x-axis: parameter values or combo index ───────────────

        function x = comboXAxis(obj)
            % COMBOXAXIS  Return x-axis values for sweep plots.
            %   If comboParamField is set, uses that parameter value.
            %   Otherwise returns 1:nCombos.
            if isempty(obj.comboParamField)
                x = 1 : obj.ana.nCombos;
                return;
            end
            x = NaN(1, obj.ana.nCombos);
            for c = 1:obj.ana.nCombos
                p = obj.ana.subject_results(c).parameters;
                if isfield(p, obj.comboParamField)
                    v = p.(obj.comboParamField);
                    if isscalar(v) && isnumeric(v)
                        x(c) = v;
                    end
                end
            end
        end

        function setComboXTicks(obj, ax, x)
            % SETCOMBOTICKS  Label x-axis ticks with parameter values.
            ax.XTick = x;
            if isempty(obj.comboParamField)
                ax.XTickLabel = arrayfun(@(c) sprintf('C%d', c), ...
                    1:obj.ana.nCombos, 'UniformOutput', false);
            else
                ax.XTickLabel = arrayfun(@(v) sprintf('%.3g', v), x, ...
                    'UniformOutput', false);
                ax.XLabel.String = obj.comboParamField;
            end
        end

        function lbl = comboShortLabel(obj, c)
            % COMBOSHORTLABEL  Short per-point annotation string.
            params = obj.ana.subject_results(c).parameters;
            flds   = fieldnames(params);
            if isempty(flds)
                lbl = sprintf('C%d', c);
                return;
            end
            parts = {};
            for i = 1:min(2, numel(flds))
                v = params.(flds{i});
                if isscalar(v) && isnumeric(v)
                    parts{end+1} = sprintf('%.4g', v); %#ok<AGROW>
                end
            end
            if isempty(parts)
                lbl = sprintf('C%d', c);
            else
                lbl = strjoin(parts, '/');
                if numel(lbl) > obj.maxLabelChars
                    lbl = lbl(1:obj.maxLabelChars);
                end
            end
        end

        % ── Generic metric puller ─────────────────────────────────

        function vals = pullMetric(obj, analyser, fieldName, seg, signalType)
            % PULLMETRIC  getSummary() dispatcher for all analysers.
            if nargin < 4, seg        = 'closed'; end
            if nargin < 5, signalType = 'diff';   end
            if isempty(seg), seg = 'closed'; end

            switch lower(analyser)
                case 'blink'
                    obj.hasBlink(true);
                    vals = obj.ana.blink.getSummary(fieldName);
                case 'signalstats'
                    obj.hasSignalStats(true);
                    vals = obj.ana.signalStats.getSummary( ...
                        fieldName, seg, signalType);
                case 'spectral'
                    obj.hasSpectral(true);
                    vals = obj.ana.spectral.getSummary(fieldName, seg);
                case 'temporal'
                    obj.hasTemporal(true);
                    vals = obj.ana.temporal.getSummary(fieldName, seg);
                otherwise
                    error('SummaryPlot:unknownAnalyser', ...
                        'Unknown analyser ''%s''. Use: blink, signalStats, spectral, temporal.', ...
                        analyser);
            end
        end

        function s = shortLabel(~, fieldName)
            % SHORTLABEL  Short axis label from dot-path field name.
            parts = strsplit(fieldName, '.');
            s     = strrep(parts{end}, '_mean', '');
            s     = strrep(s, '_', ' ');
            if numel(s) > 14, s = s(1:14); end
        end

    end % private methods

end
