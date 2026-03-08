classdef ASRPlotBase < handle
    % ASRPlotBase  Shared infrastructure for all ASR visualisation classes.
    %
    %   Parent of:  SignalPlot   — raw vs clean signals, blink events
    %               SpectralPlot — PSD, band power, topoplot metrics
    %               SummaryPlot  — parameter-sweep level comparisons
    %
    %   ── Data-access helpers (protected) ─────────────────────────────────
    %     [Xraw,Xclean] = getSegment(c, seg)
    %         Pull raw + cleaned [C×T] for combo c and segment.
    %     thresh = getThreshold(c, seg)
    %         Per-channel [C×1] blink amplitude threshold, recomputed from
    %         the method/multiplier stored in blink.results(c).summary.
    %         Returns [] if no blink analysis has been run.
    %     idx = getBlinkIdx(c, seg, role)
    %         {C×1} cell of peak sample indices.
    %         role = 'raw' | 'clean' | 'matched_raw' | 'matched_clean'
    %         Returns {} if blink analysis absent.
    %     hasBlink / hasSignalStats / hasSpectral / hasTemporal
    %         Guards so subclass methods fail cleanly with a clear message.
    %
    %   ── Figure helpers (protected) ───────────────────────────────────────
    %     fig  = newFig(name)
    %     ax   = subplotGrid(nRows, nCols, i)
    %     t    = timeAxis(T)              (0:T-1)/fs  seconds
    %     mask = windowIdx(t, t0, t1)    logical index into time vector
    %     str  = makeTitle(part1, ...)    joined with ' | '
    %     decorateAx(ax, xlab, ylab, ttl)
    %     autoSave(fig, filename)         writes to saveDir if set
    %
    %   ── Drawing primitives (protected) ──────────────────────────────────
    %     drawSignalPair(ax, t, xr, xc)
    %     drawThreshold(ax, thresh)
    %     drawBlinkMarkers(ax, tRaw, tClean, tMatched)
    %     shadeBand(ax, x, mu, sigma)
    %     addTopoColorbar(ax, label, clim)
    %
    %   ── Style properties (public, fully overridable) ─────────────────────
    %     colorRaw / colorClean / colorThresh / colorBlink
    %     colorBlinkASR / colorMatched / colorShade / colorBand
    %     lineWidth / markerSize / fontSize / figColor
    %     saveDir / saveFormat
    %
    %   ── Constructor ──────────────────────────────────────────────────────
    %     Do not instantiate ASRPlotBase directly; use a subclass:
    %       sp = SignalPlot(ana);
    %       sp = SignalPlot(ana, EEG.chanlocs);

    % ================================================================
    properties (Access = public)

        % ── Data source ────────────────────────────────────────────
        ana           % ExperimentAnalysis (already loaded)
        fs            % Hz — mirror of ana.fs
        algorithmName % string
        subjectID     % numeric

        % ── EEGLAB chanlocs for topoplot ───────────────────────────
        chanlocs = []   % [] = topo disabled

        % ── Colour palette ─────────────────────────────────────────
        colorRaw      = [0.18 0.44 0.71]  % blue        original signal
        colorClean    = [0.84 0.19 0.12]  % red         ASR-cleaned
        colorThresh   = [0.55 0.55 0.55]  % grey        threshold line
        colorBlink    = [0.10 0.10 0.10]  % near-black  raw blink peaks
        colorBlinkASR = [0.18 0.72 0.23]  % green       clean blink peaks
        colorMatched  = [0.93 0.60 0.10]  % orange      matched/persisted
        colorShade    = [0.86 0.86 0.86]  % light grey  std-band shading
        colorBand     = [0.96 0.96 0.75]  % pale yellow PSD band region

        % ── Line / marker style ────────────────────────────────────
        lineWidth  = 1.5
        markerSize = 7
        fontSize   = 11

        % ── Figure appearance ──────────────────────────────────────
        figColor = 'w'

        % ── Auto-save ──────────────────────────────────────────────
        saveDir    = ''     % '' = disabled; set to folder path to enable
        saveFormat = 'png'  % 'png' | 'pdf' | 'svg' | 'fig'

    end % properties

    % ================================================================
    methods (Access = public)

        function obj = ASRPlotBase(ana, chanlocs)
            % ASRPLOTBASE  Base constructor — called by subclass constructors.
            %   ana      : ExperimentAnalysis (already had load() called)
            %   chanlocs : (optional) EEGLAB chanlocs struct

            if nargin < 1 || isempty(ana)
                error('ASRPlotBase:noData', ...
                    'An ExperimentAnalysis instance is required.');
            end
            if isempty(ana.subject_results)
                error('ASRPlotBase:notLoaded', ...
                    'Call ana.load() before creating a plot object.');
            end

            obj.ana           = ana;
            obj.fs            = ana.fs;
            obj.algorithmName = ana.algorithmName;
            obj.subjectID     = ana.subjectID;

            if nargin >= 2 && ~isempty(chanlocs)
                obj.chanlocs = chanlocs;
            end
        end

    end % public methods

    % ================================================================
    methods (Access = protected)

        % ════════════════════════════════════════════════════════════
        %  DATA ACCESS
        % ════════════════════════════════════════════════════════════

        function [Xraw, Xclean] = getSegment(obj, c, seg)
            % GETSEGMENT  Return [C×T] raw and cleaned signal for combo c.

            obj.validateCombo(c);
            seg = obj.resolveSegment(seg);
            sr  = obj.ana.subject_results(c);

            switch seg
                case 'closed'
                    Xraw   = double(obj.ana.raw.closed);
                    Xclean = double(sr.cleanClosed);
                case 'open'
                    Xraw   = double(obj.ana.raw.open);
                    Xclean = double(sr.cleanOpen);
                case 'calibration'
                    Xraw   = double(obj.ana.raw.calibration);
                    Xclean = double(sr.cleanCalibration);
                case 'all'
                    Xraw   = double([obj.ana.raw.calibration, ...
                                     obj.ana.raw.closed, ...
                                     obj.ana.raw.open]);
                    Xclean = double([sr.cleanCalibration, ...
                                     sr.cleanClosed, ...
                                     sr.cleanOpen]);
            end

            T      = min(size(Xraw,2), size(Xclean,2));
            Xraw   = Xraw(:,   1:T);
            Xclean = Xclean(:, 1:T);
        end

        function thresh = getThreshold(obj, c, seg)
            % GETTHRESHOLD  [C×1] per-channel blink amplitude threshold.
            %   Recomputed from raw using method/multiplier stored in
            %   blink.results(c).summary.  Returns [] if unavailable.

            thresh = [];
            if ~obj.hasBlink(false) || numel(obj.ana.blink.results) < c
                return;
            end

            s   = obj.ana.blink.results(c).summary;
            k   = s.multiplier;
            seg = obj.resolveSegment(seg);

            switch s.method
                case 'mad'
                    bp = [1 10];
                    if isfield(s,'bandpass') && ~isempty(s.bandpass)
                        bp = s.bandpass;
                    end
                    [b, a] = butter(2, bp / (obj.fs/2), 'bandpass');
                    switch seg
                        case 'closed',   Xr = double(obj.ana.raw.closed);
                        case 'open',     Xr = double(obj.ana.raw.open);
                        otherwise,       Xr = double(obj.ana.raw.closed);
                    end
                    Xf     = filtfilt(b, a, Xr')';
                    thresh = k * (1.4826 * mad(Xf, 1, 2));

                otherwise  % 'mean'
                    switch seg
                        case 'closed',   Xr = double(obj.ana.raw.closed);
                        case 'open',     Xr = double(obj.ana.raw.open);
                        otherwise,       Xr = double(obj.ana.raw.closed);
                    end
                    thresh = k * mean(abs(Xr), 2);
            end
        end

        function idx = getBlinkIdx(obj, c, seg, role)
            % GETBLINKIDX  {C×1} cell of peak sample indices.
            %   role = 'raw' | 'clean' | 'matched_raw' | 'matched_clean'

            idx = {};
            if ~obj.hasBlink(false) || numel(obj.ana.blink.results) < c
                return;
            end
            seg   = obj.resolveSegment(seg);
            perCh = obj.ana.blink.results(c).perChannel;
            field = sprintf('idx_%s_%s', role, seg);
            if isfield(perCh, field)
                idx = perCh.(field);
            end
        end

        % ── Availability guards ────────────────────────────────────

        function ok = hasBlink(obj, doError)
            if nargin < 2, doError = true; end
            ok = ~isempty(obj.ana.blink) && ...
                 ~isempty(obj.ana.blink.results);
            if ~ok && doError
                error('ASRPlotBase:noBlink', ...
                    'Run ana.computeBlinks() before calling this plot method.');
            end
        end

        function ok = hasSignalStats(obj, doError)
            if nargin < 2, doError = true; end
            ok = ~isempty(obj.ana.signalStats) && ...
                 ~isempty(obj.ana.signalStats.results);
            if ~ok && doError
                error('ASRPlotBase:noSignalStats', ...
                    'Run ana.computeSignalStats() before calling this plot method.');
            end
        end

        function ok = hasSpectral(obj, doError)
            if nargin < 2, doError = true; end
            ok = isfield(obj.ana,'spectral') && ...
                 ~isempty(obj.ana.spectral) && ...
                 ~isempty(obj.ana.spectral.results);
            if ~ok && doError
                error('ASRPlotBase:noSpectral', ...
                    'Attach a computed SpectralAnalysis to ana.spectral first.');
            end
        end

        function ok = hasTemporal(obj, doError)
            if nargin < 2, doError = true; end
            ok = isfield(obj.ana,'temporal') && ...
                 ~isempty(obj.ana.temporal) && ...
                 ~isempty(obj.ana.temporal.results);
            if ~ok && doError
                error('ASRPlotBase:noTemporal', ...
                    'Attach a computed TemporalAnalysis to ana.temporal first.');
            end
        end

        % ════════════════════════════════════════════════════════════
        %  FIGURE / AXES HELPERS
        % ════════════════════════════════════════════════════════════

        function fig = newFig(obj, name)
            fig = figure('Color', obj.figColor, ...
                         'Name',  name, ...
                         'NumberTitle', 'off');
        end

        function ax = subplotGrid(~, nR, nC, i)
            ax = subplot(nR, nC, i);
        end

        function t = timeAxis(obj, T)
            t = (0 : T-1) / obj.fs;
        end

        function mask = windowIdx(~, t, t0, t1)
            mask = t >= t0 & t <= t1;
        end

        function str = makeTitle(obj, varargin)
            % MAKETITLE  'S<id>  <algo>  |  part1  |  part2  ...'
            base  = sprintf('S%d  %s', obj.subjectID, obj.algorithmName);
            extra = cellfun(@(v) num2str(v,'%g'), varargin, ...
                            'UniformOutput', false);
            str   = strjoin([{base}, extra], '  |  ');
        end

        function decorateAx(obj, ax, xLab, yLab, ttl)
            if nargin >= 3 && ~isempty(xLab), xlabel(ax, xLab); end
            if nargin >= 4 && ~isempty(yLab), ylabel(ax, yLab); end
            if nargin >= 5 && ~isempty(ttl)
                title(ax, ttl, 'Interpreter','none');
            end
            ax.FontSize = obj.fontSize;
            grid(ax, 'on');
            box(ax, 'off');
        end

        function autoSave(obj, fig, filename)
            if isempty(obj.saveDir), return; end
            if ~isfolder(obj.saveDir), mkdir(obj.saveDir); end
            fpath = fullfile(obj.saveDir, [filename, '.', obj.saveFormat]);
            saveas(fig, fpath);
            fprintf('  [saved] %s\n', fpath);
        end

        % ════════════════════════════════════════════════════════════
        %  DRAWING PRIMITIVES  (reused across subclasses)
        % ════════════════════════════════════════════════════════════

        function drawSignalPair(obj, ax, t, xr, xc)
            % DRAWSIGNALPAIR  Plot raw(blue) and clean(red) on ax.
            plot(ax, t, xr, '-', 'Color', obj.colorRaw,   ...
                 'LineWidth', obj.lineWidth, 'DisplayName', 'Raw');
            hold(ax, 'on');
            plot(ax, t, xc, '-', 'Color', obj.colorClean, ...
                 'LineWidth', obj.lineWidth, 'DisplayName', 'Clean');
        end

        function drawThreshold(obj, ax, thresh)
            % DRAWTHRESHOLD  Horizontal ±threshold dashed lines.
            if isempty(thresh), return; end
            yline(ax,  thresh, '--', 'Color', obj.colorThresh, ...
                  'LineWidth', 0.9, 'DisplayName', 'Threshold');
            if thresh > 0
                yline(ax, -thresh, '--', 'Color', obj.colorThresh, ...
                      'LineWidth', 0.9, 'HandleVisibility','off');
            end
        end

        function drawBlinkMarkers(obj, ax, tRaw, tClean, tMatched)
            % DRAWBLINKMARKERS  Scatter blink peak times at top of axis.
            %   tRaw/tClean/tMatched : time in seconds (may be empty)
            yl = ylim(ax);
            span = yl(2) - yl(1);
            if nargin >= 3 && ~isempty(tRaw)
                yv = yl(2) - 0.05*span;
                scatter(ax, tRaw, repmat(yv, size(tRaw)), ...
                    obj.markerSize*6, obj.colorBlink, '*', ...
                    'DisplayName', 'Raw blinks');
                hold(ax, 'on');
            end
            if nargin >= 4 && ~isempty(tClean)
                yv = yl(2) - 0.12*span;
                scatter(ax, tClean, repmat(yv, size(tClean)), ...
                    obj.markerSize*6, obj.colorBlinkASR, '*', ...
                    'DisplayName', 'Clean blinks');
                hold(ax, 'on');
            end
            if nargin >= 5 && ~isempty(tMatched)
                yv = yl(2) - 0.20*span;
                scatter(ax, tMatched, repmat(yv, size(tMatched)), ...
                    obj.markerSize*8, obj.colorMatched, 'o', ...
                    'LineWidth', 1.5, 'DisplayName', 'Persisted');
                hold(ax, 'on');
            end
        end

        function shadeBand(obj, ax, x, mu, sigma)
            % SHADEBAND  Translucent ±sigma fill around mu.
            px = [x(:)', fliplr(x(:)')];
            py = [mu(:)'+sigma(:)', fliplr(mu(:)'-sigma(:)')];
            fill(ax, px, py, obj.colorShade, ...
                 'EdgeColor','none', 'FaceAlpha', 0.40, ...
                 'HandleVisibility','off');
        end

        function addTopoColorbar(obj, ax, label, clim)
            % ADDTOPOCOLORBAR  Uniform colorbar styling for topoplot axes.
            cb = colorbar(ax);
            cb.Label.String   = label;
            cb.Label.FontSize = obj.fontSize;
            if nargin >= 4 && ~isempty(clim)
                caxis(ax, clim);
            end
        end

        % ════════════════════════════════════════════════════════════
        %  VALIDATION
        % ════════════════════════════════════════════════════════════

        function validateCombo(obj, c)
            if c < 1 || c > obj.ana.nCombos
                error('ASRPlotBase:badCombo', ...
                    'Combo index %d out of range [1, %d].', c, obj.ana.nCombos);
            end
        end

        function validateChannel(~, ch, C)
            if any(ch < 1) || any(ch > C)
                error('ASRPlotBase:badChannel', ...
                    'Channel index out of range [1, %d].', C);
            end
        end

        function seg = resolveSegment(~, seg)
            valid = {'closed','open','calibration','all'};
            seg   = lower(seg);
            if ~ismember(seg, valid)
                error('ASRPlotBase:badSegment', ...
                    'Segment must be one of: %s', strjoin(valid,', '));
            end
        end

        function checkTopoAvail(obj)
            if isempty(which('topoplot'))
                error('ASRPlotBase:noTopoplot', ...
                    'topoplot() not found — add EEGLAB to the path.');
            end
            if isempty(obj.chanlocs)
                error('ASRPlotBase:noChanlocs', ...
                    'chanlocs is empty — pass chanlocs to the constructor.');
            end
        end

    end % protected methods

end
