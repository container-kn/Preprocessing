classdef TemporalPlot < ASRPlotBase
    % TemporalPlot  Windowed, time-resolved metric visualisation.
    %
    %   Inherits from ASRPlotBase (colour palette, figure helpers,
    %   data-access helpers, auto-save).
    %   Reads from: ana.temporal  (TemporalAnalysis — must be computed first)
    %
    %   Complements the existing plot hierarchy:
    %     SignalPlot   — raw/clean signal traces              (time domain)
    %     SpectralPlot — PSD, band power, topoplots           (freq domain)
    %     SummaryPlot  — one scalar per combo across a sweep  (sweep level)
    %     TemporalPlot — windowed metrics over time per combo  ← this class
    %
    %   ── Plot methods ─────────────────────────────────────────────────────
    %
    %   plotSignalPanel(c, seg)
    %       3-panel: RRMSE / correlation / SNR(dB) over time for one combo.
    %       ±1 std shading across channels. Reference lines at ideal values.
    %
    %   plotSubspacePanel(c, seg)
    %       3-panel: principal angle (mean + max) / covariance Frobenius
    %       distortion / energy removed — all vs time for one combo.
    %
    %   plotEnergyCompactionPanel(c, seg)
    %       2-panel: kCompact (raw vs clean) / dimensionality reduction —
    %       shows how ASR changes the effective rank of the signal over time.
    %
    %   plotCalibrationDriftPanel(c, seg)
    %       2-panel: energy capture (raw vs clean) / angle from calibration
    %       reference — detects parameter regimes that drift away from
    %       the calibration subspace over time.
    %
    %   plotComboOverlay(fieldName, seg)
    %       All combos overlaid on one axes for a single metric, coloured
    %       by combo index. Best for identifying which combos behave differently.
    %
    %   plotComboGrid(fieldName, seg, comboIndices)
    %       Small-multiples: one subplot per combo. Each shows the same
    %       metric with its time-averaged mean as a dashed reference line.
    %       Combo title includes key parameter values.
    %
    %   plotGroupSummaryBars(seg)
    %       Grouped bar chart: one bar per combo, y = time-averaged scalar,
    %       grouped by metric group (signal / subspace / compaction / drift).
    %       Quick ranking of combos across all summary fields at once.
    %
    %   ── Properties ───────────────────────────────────────────────────────
    %     maxCombosOverlay   max combos drawn in plotComboOverlay (default 12)
    %     showMeanLine       draw dashed mean in per-combo panels (default true)
    %     showStdShade       draw ±1 std shading (default true)
    %     showRefLine        draw ideal reference line (default true)
    %     comboColormap      colormap name for multi-combo lines (default 'lines')
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       % After running analysis:
    %       ana.computeTemporal('windowSec', 2, 'hopSec', 1, 'subspaceRank', 10);
    %
    %       tp = TemporalPlot(ana);
    %       tp.saveDir = '/my/figs';          % optional auto-save
    %
    %       % Single-combo panels:
    %       tp.plotSignalPanel(1, 'closed');
    %       tp.plotSubspacePanel(1, 'closed');
    %       tp.plotEnergyCompactionPanel(1, 'closed');
    %       tp.plotCalibrationDriftPanel(1, 'closed');
    %
    %       % Cross-combo views:
    %       tp.plotComboOverlay('signal.rrmse',        'closed');
    %       tp.plotComboOverlay('subspace.angle_mean', 'closed');
    %       tp.plotComboGrid('signal.rrmse',           'closed');
    %       tp.plotGroupSummaryBars('closed');
    %
    %       % Restrict to a subset of combos:
    %       tp.plotComboGrid('signal.corr', 'closed', [1 5 10 20]);

    % ================================================================
    properties (Access = public)
        maxCombosOverlay = 12     % cap for plotComboOverlay (avoids overplotting)
        showMeanLine     = true   % dashed mean reference in per-combo panels
        showStdShade     = true   % ±1 std shading across channels
        showRefLine      = true   % horizontal ideal value reference (0, 1, etc.)
        comboColormap    = 'lines' % colormap name for multi-combo overlays
    end

    % ================================================================
    methods (Access = public)

        % ── Constructor ──────────────────────────────────────────────
        function obj = TemporalPlot(ana, chanlocs)
            % TEMPORALPLOT
            %   ana      : ExperimentAnalysis (ana.computeTemporal() must
            %              have been called before any plot method)
            %   chanlocs : (optional) EEGLAB chanlocs — not used by this
            %              class but passed through to ASRPlotBase for
            %              consistency with the rest of the hierarchy
            if nargin < 2, chanlocs = []; end
            obj = obj@ASRPlotBase(ana, chanlocs);
        end

        % ════════════════════════════════════════════════════════════
        %  1. SIGNAL PANEL — RRMSE / corr / SNR over time
        % ════════════════════════════════════════════════════════════
        function fig = plotSignalPanel(obj, c, seg)
            % PLOTSIGNALPANEL  3-panel windowed signal quality for one combo.
            %
            %   Panels: RRMSE (lower = better, ref=0),
            %           Pearson correlation (higher = better, ref=1),
            %           SNR dB (higher = better, ref=NaN).
            %
            %   c   : combo index  (default 1)
            %   seg : 'closed' | 'open'  (default 'closed')
            %
            %   Example:
            %     tp.plotSignalPanel(1, 'closed');

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end
            obj.hasTemporal(true);
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            r    = obj.ana.temporal.results(c).segment.(seg);
            t    = r.time;

            metrics = {
                'rrmse',  'RRMSE',        0,   'lower = better';
                'corr',   'Correlation',  1,   'higher = better';
                'snrDB',  'SNR (dB)',     NaN, 'higher = better'};

            fig = obj.newFig(obj.makeTitle( ...
                'Signal panel', seg, sprintf('combo %d', c)));

            for pi = 1:3
                fld  = metrics{pi, 1};
                ylab = metrics{pi, 2};
                ref  = metrics{pi, 3};

                v = r.signal.(fld);

                ax = obj.subplotGrid(1, 3, pi);
                obj.drawTemporalTrace(ax, t, v, ref, ylab);
                obj.decorateAx(ax, 'Time (s)', ylab, ...
                    sprintf('%s  |  %s  |  %s', ylab, seg, metrics{pi,4}));
            end

            sgtitle(obj.makeTitle('Signal quality over time', seg, ...
                sprintf('combo %d', c)), ...
                'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('temporal_signal_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  2. SUBSPACE PANEL — angles / covariance / energy removed
        % ════════════════════════════════════════════════════════════
        function fig = plotSubspacePanel(obj, c, seg)
            % PLOTSUBSPACEPANEL  3-panel windowed subspace geometry for one combo.
            %
            %   Panels: mean principal angle (raw vs clean subspace),
            %           max principal angle,
            %           covariance Frobenius distortion (norm).
            %
            %   c   : combo index  (default 1)
            %   seg : 'closed' | 'open'  (default 'closed')
            %
            %   Example:
            %     tp.plotSubspacePanel(1, 'closed');

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end
            obj.hasTemporal(true);
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            r = obj.ana.temporal.results(c).segment.(seg);
            t = r.time;

            metrics = {
                'angle_mean',      'Mean principal angle (°)',  0;
                'angle_max',       'Max principal angle (°)',   0;
                'cov_fro_norm',    'Cov Frobenius distortion',  0};

            fig = obj.newFig(obj.makeTitle( ...
                'Subspace panel', seg, sprintf('combo %d', c)));

            for pi = 1:3
                fld  = metrics{pi, 1};
                ylab = metrics{pi, 2};
                ref  = metrics{pi, 3};

                if ~isfield(r.subspace, fld), continue; end
                v = r.subspace.(fld);

                ax = obj.subplotGrid(1, 3, pi);
                obj.drawTemporalTrace(ax, t, v, ref, ylab);
                obj.decorateAx(ax, 'Time (s)', ylab, ...
                    sprintf('%s  |  %s', ylab, seg));
            end

            sgtitle(obj.makeTitle('Subspace geometry over time', seg, ...
                sprintf('combo %d', c)), ...
                'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('temporal_subspace_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  3. ENERGY COMPACTION PANEL — kCompact / dim reduction
        % ════════════════════════════════════════════════════════════
        function fig = plotEnergyCompactionPanel(obj, c, seg)
            % PLOTENERGYCOMPACTIONPANEL  2-panel energy compaction over time.
            %
            %   Panel 1: kCompact raw (blue) vs clean (red) — number of
            %            dimensions needed to capture thresh% of energy.
            %   Panel 2: dimensionality reduction (raw - clean kCompact).
            %            Positive = ASR reduced rank; negative = increased.
            %
            %   c   : combo index  (default 1)
            %   seg : 'closed' | 'open'  (default 'closed')
            %
            %   Example:
            %     tp.plotEnergyCompactionPanel(1, 'closed');

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end
            obj.hasTemporal(true);
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            r   = obj.ana.temporal.results(c).segment.(seg);
            t   = r.time;
            ecp = r.energyCompaction;

            thresh = obj.ana.temporal.compactionThresh * 100;

            fig = obj.newFig(obj.makeTitle( ...
                'Energy compaction', seg, sprintf('combo %d', c)));

            % Panel 1: raw vs clean kCompact overlaid
            ax1 = obj.subplotGrid(1, 2, 1);
            hold(ax1, 'on');
            plot(ax1, t, ecp.kCompact_raw,   '-', ...
                 'Color', obj.colorRaw,   'LineWidth', obj.lineWidth, ...
                 'DisplayName', 'Raw');
            plot(ax1, t, ecp.kCompact_clean, '-', ...
                 'Color', obj.colorClean, 'LineWidth', obj.lineWidth, ...
                 'DisplayName', 'Clean');
            legend(ax1, 'Location','best');
            if obj.showMeanLine
                yline(ax1, mean(ecp.kCompact_raw,   'omitnan'), '--', ...
                      'Color', obj.colorRaw,   'LineWidth', 0.8, ...
                      'HandleVisibility','off');
                yline(ax1, mean(ecp.kCompact_clean, 'omitnan'), '--', ...
                      'Color', obj.colorClean, 'LineWidth', 0.8, ...
                      'HandleVisibility','off');
            end
            obj.decorateAx(ax1, 'Time (s)', ...
                sprintf('k for %.0f%% energy', thresh), ...
                sprintf('kCompact  |  %s  |  combo %d', seg, c));

            % Panel 2: dimensionality reduction
            ax2 = obj.subplotGrid(1, 2, 2);
            dimRed = ecp.kCompact_raw - ecp.kCompact_clean;
            obj.drawTemporalTrace(ax2, t, dimRed, 0, 'Dim reduction');
            obj.decorateAx(ax2, 'Time (s)', '\Deltadims (raw - clean)', ...
                sprintf('Dimensionality reduction  |  %s  |  combo %d', seg, c));

            sgtitle(obj.makeTitle('Energy compaction over time', seg, ...
                sprintf('combo %d', c)), ...
                'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('temporal_compaction_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  4. CALIBRATION DRIFT PANEL — energy capture / angle from ref
        % ════════════════════════════════════════════════════════════
        function fig = plotCalibrationDriftPanel(obj, c, seg)
            % PLOTCALIBRATIONDRIFTPANEL  2-panel calibration drift over time.
            %
            %   Panel 1: energy capture (fraction) — raw (blue) vs clean (red).
            %            How much of the window's energy is in the calibration
            %            subspace. ASR should push clean closer to calibration.
            %   Panel 2: angle from calibration reference (°) — raw vs clean.
            %            Lower clean angle = better alignment to calibration.
            %
            %   Requires calibration segment to have been present at compute time.
            %
            %   c   : combo index  (default 1)
            %   seg : 'closed' | 'open'  (default 'closed')
            %
            %   Example:
            %     tp.plotCalibrationDriftPanel(1, 'closed');

            if nargin < 2, c   = 1;        end
            if nargin < 3, seg = 'closed'; end
            obj.hasTemporal(true);
            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);

            r   = obj.ana.temporal.results(c).segment.(seg);
            t   = r.time;
            cal = r.calibrationDrift;

            if ~isfield(cal, 'energyCapture_raw') || isempty(cal.energyCapture_raw)
                warning('TemporalPlot:noCalibDrift', ...
                    ['calibrationDrift fields empty for S%d combo %d.\n' ...
                     'Requires a calibration segment at TemporalAnalysis.compute() time.'], ...
                    obj.subjectID, c);
                fig = [];
                return;
            end

            fig = obj.newFig(obj.makeTitle( ...
                'Calibration drift', seg, sprintf('combo %d', c)));

            % Panel 1: energy capture
            ax1 = obj.subplotGrid(1, 2, 1);
            hold(ax1, 'on');
            plot(ax1, t, cal.energyCapture_raw,   '-', ...
                 'Color', obj.colorRaw,   'LineWidth', obj.lineWidth, ...
                 'DisplayName', 'Raw');
            plot(ax1, t, cal.energyCapture_clean, '-', ...
                 'Color', obj.colorClean, 'LineWidth', obj.lineWidth, ...
                 'DisplayName', 'Clean');
            if obj.showRefLine
                yline(ax1, 1, ':', 'Color', obj.colorThresh, ...
                      'LineWidth', 0.8, 'Label','ideal=1', ...
                      'LabelHorizontalAlignment','left', ...
                      'HandleVisibility','off');
            end
            legend(ax1, 'Location','best');
            obj.decorateAx(ax1, 'Time (s)', 'Energy capture (fraction)', ...
                sprintf('Energy in calibration subspace  |  %s  |  combo %d', seg, c));

            % Panel 2: angle from calibration reference
            ax2 = obj.subplotGrid(1, 2, 2);
            hold(ax2, 'on');
            plot(ax2, t, cal.angleFromRef_raw,   '-', ...
                 'Color', obj.colorRaw,   'LineWidth', obj.lineWidth, ...
                 'DisplayName', 'Raw');
            plot(ax2, t, cal.angleFromRef_clean, '-', ...
                 'Color', obj.colorClean, 'LineWidth', obj.lineWidth, ...
                 'DisplayName', 'Clean');
            if obj.showRefLine
                yline(ax2, 0, ':', 'Color', obj.colorThresh, ...
                      'LineWidth', 0.8, 'Label','ideal=0', ...
                      'LabelHorizontalAlignment','left', ...
                      'HandleVisibility','off');
            end
            legend(ax2, 'Location','best');
            obj.decorateAx(ax2, 'Time (s)', 'Angle from ref (°)', ...
                sprintf('Subspace angle from calibration  |  %s  |  combo %d', seg, c));

            sgtitle(obj.makeTitle('Calibration drift over time', seg, ...
                sprintf('combo %d', c)), ...
                'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('temporal_caldrift_c%d_%s', c, seg));
        end

        % ════════════════════════════════════════════════════════════
        %  5. COMBO OVERLAY — all combos, one metric
        % ════════════════════════════════════════════════════════════
        function fig = plotComboOverlay(obj, fieldName, seg, comboIndices)
            % PLOTCOMBOOVERLAY  All combos overlaid for one windowed metric.
            %
            %   Each combo is one coloured line. Useful for immediately seeing
            %   which parameter settings produce outlier behaviour over time.
            %
            %   fieldName    : 'group.field'  e.g. 'signal.rrmse',
            %                  'subspace.angle_mean', 'energyCompaction.kCompact_raw'
            %   seg          : 'closed' | 'open'  (default 'closed')
            %   comboIndices : subset of combos to overlay (default: all, capped
            %                  at maxCombosOverlay)
            %
            %   Example:
            %     tp.plotComboOverlay('signal.rrmse',        'closed');
            %     tp.plotComboOverlay('subspace.angle_mean', 'open', [1 5 10]);

            if nargin < 3, seg          = 'closed';            end
            if nargin < 4, comboIndices = 1:obj.ana.nCombos;   end
            obj.hasTemporal(true);
            seg = obj.resolveSegment(seg);

            % Cap to maxCombosOverlay
            if numel(comboIndices) > obj.maxCombosOverlay
                warning('TemporalPlot:tooManyCombos', ...
                    'Showing first %d of %d combos. Set tp.maxCombosOverlay to change.', ...
                    obj.maxCombosOverlay, numel(comboIndices));
                comboIndices = comboIndices(1:obj.maxCombosOverlay);
            end

            nDraw = numel(comboIndices);
            cmap  = colormap_n(obj.comboColormap, nDraw);

            fig = obj.newFig(obj.makeTitle(strrep(fieldName,'.','  '), seg, 'all combos'));
            ax  = axes(fig);
            hold(ax, 'on');

            legLabels = cell(nDraw, 1);

            for di = 1:nDraw
                c   = comboIndices(di);
                [t, v] = obj.getField(c, fieldName, seg);
                if isempty(t), continue; end

                plot(ax, t, v, '-', ...
                     'Color', cmap(di,:), ...
                     'LineWidth', obj.lineWidth * 0.9, ...
                     'DisplayName', obj.comboShortLabel(c));
                legLabels{di} = obj.comboShortLabel(c);
            end

            legend(ax, 'Location','northeastoutside', ...
                   'FontSize', obj.fontSize-2, 'Interpreter','none');
            ylab = strrep(fieldName, '.', ' — ');
            obj.decorateAx(ax, 'Time (s)', ylab, ...
                obj.makeTitle(ylab, seg, sprintf('%d combos', nDraw)));
            obj.autoSave(fig, sprintf('overlay_%s_%s', ...
                strrep(fieldName,'.','_'), seg));
        end

        % ════════════════════════════════════════════════════════════
        %  6. COMBO GRID — small-multiples, one panel per combo
        % ════════════════════════════════════════════════════════════
        function fig = plotComboGrid(obj, fieldName, seg, comboIndices)
            % PLOTCOMBOGRID  One subplot per combo for one windowed metric.
            %
            %   Each subplot shows the time series for that combo with:
            %     - the trace coloured by combo index
            %     - a dashed line at the time-averaged mean (if showMeanLine)
            %     - an ideal reference line (if showRefLine)
            %   Title: short parameter label for the combo.
            %
            %   fieldName    : 'group.field'  e.g. 'signal.rrmse'
            %   seg          : 'closed' | 'open'  (default 'closed')
            %   comboIndices : subset to show (default: all)
            %
            %   Example:
            %     tp.plotComboGrid('signal.rrmse', 'closed');
            %     tp.plotComboGrid('signal.corr',  'closed', [1 5 10 20]);

            if nargin < 3, seg          = 'closed';            end
            if nargin < 4, comboIndices = 1:obj.ana.nCombos;   end
            obj.hasTemporal(true);
            seg = obj.resolveSegment(seg);

            nDraw = numel(comboIndices);
            nC    = ceil(sqrt(nDraw));
            nR    = ceil(nDraw / nC);
            cmap  = colormap_n(obj.comboColormap, nDraw);
            ylab  = strrep(fieldName, '.', ' — ');

            fig = obj.newFig(obj.makeTitle(ylab, seg, 'combo grid'));

            for di = 1:nDraw
                c      = comboIndices(di);
                [t, v] = obj.getField(c, fieldName, seg);
                if isempty(t), continue; end

                ax = obj.subplotGrid(nR, nC, di);
                hold(ax, 'on');

                plot(ax, t, v, '-', ...
                     'Color', cmap(di,:), 'LineWidth', obj.lineWidth * 0.9, ...
                     'HandleVisibility','off');

                if obj.showMeanLine && ~all(isnan(v))
                    yline(ax, mean(v,'omitnan'), '--', ...
                          'Color', [0.3 0.3 0.3], 'LineWidth', 0.8, ...
                          'HandleVisibility','off');
                end

                ax.FontSize = obj.fontSize - 2;
                title(ax, obj.comboShortLabel(c), ...
                      'FontSize', obj.fontSize-2, 'Interpreter','none');
                xlabel(ax, 'Time (s)');
                ylabel(ax, ylab);
                grid(ax, 'on'); box(ax, 'off');
            end

            sgtitle(obj.makeTitle(ylab, seg, 'per combo'), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('grid_%s_%s', ...
                strrep(fieldName,'.','_'), seg));
        end

        % ════════════════════════════════════════════════════════════
        %  7. GROUP SUMMARY BARS — all combos × all summary scalars
        % ════════════════════════════════════════════════════════════
        function fig = plotGroupSummaryBars(obj, seg)
            % PLOTGROUPSUMMARYBARS  Bar chart across combos for key summary scalars.
            %
            %   One figure with four panels — one per metric group:
            %     Signal   : rrmse_mean, corr_mean
            %     Subspace : angle_mean_mean, cov_fro_norm_mean, energy_removed_mean
            %     Compaction: dimReduction_mean, kCompact_clean_mean
            %     CalibDrift: energyCaptureGain_mean, angleFromRefGain_mean
            %
            %   Each panel: combos on x-axis, metric value on y-axis.
            %   Multiple metrics in a panel are drawn as grouped bars.
            %
            %   seg : 'closed' | 'open'  (default 'closed')
            %
            %   Example:
            %     tp.plotGroupSummaryBars('closed');

            if nargin < 2, seg = 'closed'; end
            obj.hasTemporal(true);
            seg = obj.resolveSegment(seg);

            nC = obj.ana.nCombos;
            x  = 1:nC;

            groups = { ...
                'Signal quality', ...
                    {{'signal','rrmse_mean'},   'RRMSE'; ...
                     {'signal','corr_mean'},    'Corr'};
                'Subspace geometry', ...
                    {{'subspace','angle_mean_mean'},     'Angle (mean)'; ...
                     {'subspace','cov_fro_norm_mean'},   'Cov Fro'; ...
                     {'subspace','energy_removed_mean'}, 'E removed'};
                'Energy compaction', ...
                    {{'energyCompaction','kCompact_clean_mean'}, 'kClean'; ...
                     {'energyCompaction','dimReduction_mean'},   'DimRed'};
                'Calibration drift', ...
                    {{'calibrationDrift','energyCaptureGain_mean'},  'ECapGain'; ...
                     {'calibrationDrift','angleFromRefGain_mean'},   'AngleGain'} ...
            };

            nGroups = size(groups, 1);
            fig = obj.newFig(obj.makeTitle('Group summary bars', seg));

            for gi = 1:nGroups
                groupTitle   = groups{gi, 1};
                metricSpecs  = groups{gi, 2};
                nM           = size(metricSpecs, 1);

                data = NaN(nC, nM);
                for mi = 1:nM
                    grp = metricSpecs{mi,1}{1};
                    fld = metricSpecs{mi,1}{2};
                    for c = 1:nC
                        try
                            s = obj.ana.temporal.results(c).segment.(seg).summary;
                            data(c, mi) = s.(grp).(fld);
                        catch; end
                    end
                end

                ax = obj.subplotGrid(2, 2, gi);
                hold(ax, 'on');

                if nM == 1
                    bar(ax, x, data(:,1), 'FaceColor', obj.colorRaw);
                else
                    b = bar(ax, x, data, 'grouped');
                    colours = lines(nM);
                    for mi = 1:nM
                        b(mi).FaceColor = colours(mi,:);
                    end
                    legLabels = metricSpecs(:, 2);
                    legend(ax, legLabels{:}, 'Location','northeastoutside', ...
                           'FontSize', obj.fontSize-3);
                end

                obj.decorateAx(ax, 'Combo', 'Value', ...
                    sprintf('%s  |  %s', groupTitle, seg));
                ax.XTick     = x;
                ax.XTickLabel = arrayfun(@(c) sprintf('%d',c), x, 'UniformOutput',false);
                ax.FontSize  = obj.fontSize - 1;
            end

            sgtitle(obj.makeTitle('Temporal summary — all groups', seg), ...
                    'FontSize', obj.fontSize+1, 'Interpreter','none');
            obj.autoSave(fig, sprintf('temporal_summary_bars_%s', seg));
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        % ── drawTemporalTrace ─────────────────────────────────────────
        function drawTemporalTrace(obj, ax, t, v, refVal, ylab)
            % DRAWTEMPORALTRACE  Single-series trace with optional mean + ref.
            %   Reused by plotSignalPanel, plotSubspacePanel, plotComboGrid.

            hold(ax, 'on');

            % Std shading — not available from summary scalar, skip
            % (std across channels is not stored in time-series; shading
            %  is reserved for plotComboOverlay where we have many combos)

            % Main trace
            plot(ax, t, v, '-', ...
                 'Color',     obj.colorRaw, ...
                 'LineWidth', obj.lineWidth, ...
                 'HandleVisibility','off');

            % Mean dashed line
            if obj.showMeanLine && ~all(isnan(v))
                mu = mean(v, 'omitnan');
                yline(ax, mu, '--', 'Color', [0.35 0.35 0.35], ...
                      'LineWidth', 0.9, 'Label', sprintf('mean=%.3g', mu), ...
                      'LabelHorizontalAlignment','left', ...
                      'HandleVisibility','off');
            end

            % Reference line
            if obj.showRefLine && ~isnan(refVal)
                yline(ax, refVal, ':', 'Color', obj.colorThresh, ...
                      'LineWidth', 0.9, 'Label', sprintf('ideal=%.3g', refVal), ...
                      'LabelHorizontalAlignment','right', ...
                      'HandleVisibility','off');
            end

            grid(ax, 'on'); box(ax, 'off');
            ax.FontSize = obj.fontSize;
        end

        % ── getField ─────────────────────────────────────────────────
        function [t, v] = getField(obj, c, fieldName, seg)
            % GETFIELD  Extract time axis + one series for combo c.
            %   fieldName: 'group.field'  e.g. 'signal.rrmse'

            t = []; v = [];
            try
                [t, v] = obj.ana.temporal.getTimeSeries(c, fieldName, seg);
            catch e
                warning('TemporalPlot:missingField', ...
                    'S%d combo %d  fieldName="%s"  seg="%s": %s', ...
                    obj.subjectID, c, fieldName, seg, e.message);
            end
        end

        % ── comboShortLabel ──────────────────────────────────────────
        function lbl = comboShortLabel(obj, c)
            % COMBOSHORTLABEL  Short key=value summary of a combo's params.
            %   Shows up to 2 most-varied scalar parameters.

            params = obj.ana.subject_results(c).parameters;
            if isempty(params)
                lbl = sprintf('c%d', c);
                return;
            end
            fields = fieldnames(params);
            parts  = cell(1, numel(fields));
            for f = 1:numel(fields)
                v = params.(fields{f});
                if isscalar(v) && isnumeric(v)
                    parts{f} = sprintf('%s=%.4g', fields{f}, v);
                elseif ischar(v)
                    parts{f} = sprintf('%s=%s', fields{f}, v);
                end
            end
            parts(cellfun(@isempty, parts)) = [];
            if isempty(parts)
                lbl = sprintf('c%d', c);
            else
                lbl = strjoin(parts(1:min(2,end)), ' ');
            end
        end

    end % private methods

end % classdef

% ── module-level helper: N colours from a named colormap ─────────────────
function cmap = colormap_n(name, n)
    cmap = feval(name, max(n, 1));
    if size(cmap, 1) < n
        cmap = repmat(cmap, ceil(n/size(cmap,1)), 1);
    end
    cmap = cmap(1:n, :);
end
