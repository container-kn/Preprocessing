classdef BlinkAnalysis < handle
    % BlinkAnalysis - Eye blink detection, reduction, and persistence metrics.
    %
    %   Self-contained analyser — receives raw and cleaned signals, runs
    %   detection across all parameter combinations, and exposes results.
    %
    %   ── Detection methods ───────────────────────────────────────────────
    %     'mean' (default) — threshold = k * mean(|raw|)  per channel
    %     'mad'            — bandpass [1 10] Hz, threshold = k * 1.4826*MAD(raw)
    %   Threshold always fixed from RAW only — identical across all combos.
    %
    %   ── Asymmetric peak separation ──────────────────────────────────────
    %   Raw and clean signals benefit from different minimum peak distances:
    %     minDistRaw   (default 0.2 s) — raw has sharp, well-separated peaks
    %     minDistClean (default 0.5 s) — cleaned signal may have ringing near
    %                                    blinks; stricter separation prevents
    %                                    one blink from generating two clean peaks
    %
    %   ── Persistence — bipartite matching ────────────────────────────────
    %   A blink is "persistent" if it survived ASR: a clean peak can be
    %   matched to an original raw peak within ±toleranceSec.
    %
    %   WHY NAIVE MATCHING FAILS:
    %   The naive check "for each clean blink, is there ANY raw blink within tol?"
    %   lets multiple clean blinks claim the same raw blink, making the
    %   persistent count exceed the raw count — a nonsensical result.
    %
    %   CORRECT APPROACH — greedy 1-D bipartite matching:
    %     1. Shift clean indices to raw time base (+offset, default 0)
    %     2. For each clean blink (sorted in time), find the nearest
    %        UNCLAIMED raw blink within ±toleranceSamples
    %     3. Claim that raw blink — no other clean blink can take it
    %     4. nMatched = total successful claims
    %   This guarantees nMatched <= min(n_raw, n_clean) and is optimal
    %   for 1-D sorted events with a symmetric tolerance window.
    %
    %   Offset: our ASR output is already time-aligned with raw (the P+1:end
    %   carry-stripping ensures this), so offset = 0 by default. Expose it
    %   for pipelines that strip lookahead externally.
    %
    %   ── Per-channel results (results(k).perChannel) ─────────────────────
    %     nBlinks_raw/clean_closed/open       [C x 1]
    %     nMatched_closed/open                [C x 1]  bipartite match count
    %     nRemoved_closed/open                [C x 1]  raw - matched
    %     nNewArtefacts_closed/open           [C x 1]  clean - matched
    %     idx_raw/clean_closed/open           {C x 1}  peak sample indices
    %     idx_matched_raw/clean_closed/open   {C x 1}  matched peak indices
    %
    %   ── Summary scalars (results(k).summary) ────────────────────────────
    %     totalBlinks_raw/clean_closed/open
    %     meanBlinks_raw/clean_closed/open
    %     blinkReductionRatio_closed/open     1 - clean/raw, clamped [0,1]
    %     blinkReductionPct_closed/open       ratio * 100
    %     blinkReductionAbs_closed/open       raw - clean  (integer)
    %     totalMatched_closed/open            sum of nMatched across channels
    %     persistenceRate_closed/open         nMatched/n_raw  (fraction NOT removed)
    %     removalRate_closed/open             1 - persistenceRate
    %     artefactRate_closed/open            nNewArtefacts/n_clean
    %     valid_closed/open                   physiological plausibility flag
    %     amplified_closed/open               clean > raw flag
    %
    %   ── Usage ────────────────────────────────────────────────────────────
    %       bl = BlinkAnalysis(raw, subject_results, fs, subjectID, algo);
    %       bl.compute('method','mean','toleranceSec',0.02);
    %       bl.getSummary('persistenceRate_closed')
    %       bl.getSummary('removalRate_closed')
    %       bl.getSummary('blinkReductionRatio_closed')
    %       bl.paramTable()

    % ================================================================
    properties (Access = public)
        raw
        subject_results
        fs
        subjectID
        algorithmName
        nCombos
        results
    end

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = BlinkAnalysis(raw, subject_results, fs, subjectID, algorithmName)
            if nargin < 4, subjectID     = 0;         end
            if nargin < 5, algorithmName = 'unknown'; end

            obj.raw             = raw;
            obj.subject_results = subject_results;
            obj.fs              = fs;
            obj.subjectID       = subjectID;
            obj.algorithmName   = algorithmName;
            obj.nCombos         = numel(subject_results);
            obj.results         = struct();
        end

        % ── compute ──────────────────────────────────────────────────
        function compute(obj, varargin)
            % COMPUTE - Detect blinks and persistence across all combos.
            %
            %   Name-value options:
            %     'method'        'mean'   Detection: 'mean' | 'mad'
            %     'bandpass'      [1 10]   Hz band for 'mad' method only
            %     'multiplier'    []       Threshold k  (default: 6 mean / 5 mad)
            %     'minDistRaw'    []       Min sec between raw peaks  (default 0.2)
            %     'minDistClean'  []       Min sec between clean peaks (default 0.5)
            %     'maxBlinkRate'  0.5      Max plausible blinks/sec (sanity check)
            %     'toleranceSec'  0.02     Persistence match window ±sec
            %     'offset'        0        Samples added to clean indices before
            %                              matching (0 for our pipeline; non-zero
            %                              only if clean was trimmed externally)

            p = inputParser;
            p.KeepUnmatched = false;
            addParameter(p, 'method',       'mean', @(x) ismember(x, {'mad','mean'}));
            addParameter(p, 'bandpass',     [1 10], @(x) isnumeric(x) && numel(x)==2);
            addParameter(p, 'multiplier',   [],     @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'minDistRaw',   [],     @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'minDistClean', [],     @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'maxBlinkRate', 0.5,    @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'toleranceSec', 0.02,   @(x) isnumeric(x) && isscalar(x) && x >= 0);
            addParameter(p, 'offset',       0,      @(x) isnumeric(x) && isscalar(x));
            parse(p, varargin{:});

            method       = p.Results.method;
            bp           = p.Results.bandpass;
            maxBlinkRate = p.Results.maxBlinkRate;
            toleranceSec = p.Results.toleranceSec;
            offset       = round(p.Results.offset);

            % Per-method multiplier defaults
            if isempty(p.Results.multiplier)
                k = struct('mad', 5, 'mean', 6).(method);
            else
                k = p.Results.multiplier;
            end

            % Asymmetric minPeakDistance defaults
            if isempty(p.Results.minDistRaw)
                minDistRaw = struct('mad', 0.5, 'mean', 0.2).(method);
            else
                minDistRaw = p.Results.minDistRaw;
            end
            if isempty(p.Results.minDistClean)
                minDistClean = 0.5;
            else
                minDistClean = p.Results.minDistClean;
            end

            fprintf('BlinkAnalysis [%s]: S%d %s (%d combo(s))...\n', ...
                method, obj.subjectID, obj.algorithmName, obj.nCombos);

            % Pre-compute thresholds from RAW once
            switch method
                case 'mad'
                    [b, a]       = butter(2, bp / (obj.fs / 2), 'bandpass');
                    rawCl_filt   = filtfilt(b, a, double(obj.raw.closed)')';
                    rawOp_filt   = filtfilt(b, a, double(obj.raw.open)')';
                    threshClosed = k * (1.4826 * mad(rawCl_filt, 1, 2));
                    threshOpen   = k * (1.4826 * mad(rawOp_filt, 1, 2));
                case 'mean'
                    b = []; a = [];
                    threshClosed = k * mean(abs(double(obj.raw.closed)), 2);
                    threshOpen   = k * mean(abs(double(obj.raw.open)),   2);
            end

            durClosed  = size(obj.raw.closed, 2) / obj.fs;
            durOpen    = size(obj.raw.open,   2) / obj.fs;
            maxClosed  = maxBlinkRate * durClosed;
            maxOpen    = maxBlinkRate * durOpen;
            tolSamples = round(toleranceSec * obj.fs);

            for c = 1:obj.nCombos

                cleanClosed = obj.subject_results(c).cleanClosed;
                cleanOpen   = obj.subject_results(c).cleanOpen;

                % --- Peak detection (asymmetric minDist) ---
                [n_raw_cl, n_cln_cl, idx_raw_cl, idx_cln_cl] = ...
                    BlinkAnalysis.detectBlinks( ...
                        obj.raw.closed, cleanClosed, threshClosed, ...
                        minDistRaw, minDistClean, obj.fs, b, a, method);

                [n_raw_op, n_cln_op, idx_raw_op, idx_cln_op] = ...
                    BlinkAnalysis.detectBlinks( ...
                        obj.raw.open, cleanOpen, threshOpen, ...
                        minDistRaw, minDistClean, obj.fs, b, a, method);

                % --- Bipartite persistence matching ---
                [nMatch_cl, idxMR_cl, idxMC_cl] = BlinkAnalysis.matchBlinks( ...
                    idx_raw_cl, idx_cln_cl, tolSamples, offset);

                [nMatch_op, idxMR_op, idxMC_op] = BlinkAnalysis.matchBlinks( ...
                    idx_raw_op, idx_cln_op, tolSamples, offset);

                % --- Sanity warnings ---
                warn_cl = find((n_raw_cl > maxClosed) | (n_cln_cl > maxClosed));
                warn_op = find((n_raw_op > maxOpen)   | (n_cln_op > maxOpen));
                if ~isempty(warn_cl)
                    fprintf(['  [BlinkAnalysis] WARNING: Closed-eyes channels %s ' ...
                        'exceed max plausible blink count (%.0f).\n'], ...
                        mat2str(warn_cl), maxClosed);
                end
                if ~isempty(warn_op)
                    fprintf(['  [BlinkAnalysis] WARNING: Open-eyes channels %s ' ...
                        'exceed max plausible blink count (%.0f).\n'], ...
                        mat2str(warn_op), maxOpen);
                end

                % --- Per-channel struct ---
                perCh = struct( ...
                    'nBlinks_raw_closed',        n_raw_cl,      ...
                    'nBlinks_clean_closed',       n_cln_cl,      ...
                    'nMatched_closed',            nMatch_cl,     ...
                    'nRemoved_closed',            n_raw_cl - nMatch_cl, ...
                    'nNewArtefacts_closed',       n_cln_cl - nMatch_cl, ...
                    'idx_raw_closed',             {idx_raw_cl},  ...
                    'idx_clean_closed',           {idx_cln_cl},  ...
                    'idx_matched_raw_closed',     {idxMR_cl},    ...
                    'idx_matched_clean_closed',   {idxMC_cl},    ...
                    'nBlinks_raw_open',           n_raw_op,      ...
                    'nBlinks_clean_open',         n_cln_op,      ...
                    'nMatched_open',              nMatch_op,     ...
                    'nRemoved_open',              n_raw_op - nMatch_op, ...
                    'nNewArtefacts_open',         n_cln_op - nMatch_op, ...
                    'idx_raw_open',               {idx_raw_op},  ...
                    'idx_clean_open',             {idx_cln_op},  ...
                    'idx_matched_raw_open',       {idxMR_op},    ...
                    'idx_matched_clean_open',     {idxMC_op});

                % --- Summary struct ---
                summary = struct( ...
                    'totalBlinks_raw_closed',     sum(n_raw_cl),  ...
                    'totalBlinks_clean_closed',   sum(n_cln_cl),  ...
                    'meanBlinks_raw_closed',      mean(n_raw_cl), ...
                    'meanBlinks_clean_closed',    mean(n_cln_cl), ...
                    'totalBlinks_raw_open',       sum(n_raw_op),  ...
                    'totalBlinks_clean_open',     sum(n_cln_op),  ...
                    'meanBlinks_raw_open',        mean(n_raw_op), ...
                    'meanBlinks_clean_open',      mean(n_cln_op), ...
                    'blinkReductionRatio_closed', BlinkAnalysis.reductionRatio(n_raw_cl, n_cln_cl), ...
                    'blinkReductionRatio_open',   BlinkAnalysis.reductionRatio(n_raw_op, n_cln_op), ...
                    'blinkReductionPct_closed',   100*BlinkAnalysis.reductionRatio(n_raw_cl, n_cln_cl), ...
                    'blinkReductionPct_open',     100*BlinkAnalysis.reductionRatio(n_raw_op, n_cln_op), ...
                    'blinkReductionAbs_closed',   sum(n_raw_cl) - sum(n_cln_cl), ...
                    'blinkReductionAbs_open',     sum(n_raw_op) - sum(n_cln_op), ...
                    'totalMatched_closed',        sum(nMatch_cl), ...
                    'totalMatched_open',          sum(nMatch_op), ...
                    'persistenceRate_closed',     BlinkAnalysis.persistenceRate(n_raw_cl, nMatch_cl), ...
                    'persistenceRate_open',       BlinkAnalysis.persistenceRate(n_raw_op, nMatch_op), ...
                    'removalRate_closed',         BlinkAnalysis.removalRate(n_raw_cl, nMatch_cl),    ...
                    'removalRate_open',           BlinkAnalysis.removalRate(n_raw_op, nMatch_op),    ...
                    'artefactRate_closed',        BlinkAnalysis.artefactRate(n_cln_cl, nMatch_cl),   ...
                    'artefactRate_open',          BlinkAnalysis.artefactRate(n_cln_op, nMatch_op),   ...
                    'parameters',    obj.subject_results(c).parameters, ...
                    'method',        method,       ...
                    'multiplier',    k,            ...
                    'minDistRaw',    minDistRaw,   ...
                    'minDistClean',  minDistClean, ...
                    'toleranceSec',  toleranceSec, ...
                    'offset',        offset);

                if strcmp(method, 'mad')
                    summary.bandpass = bp;
                end

                summary.valid_closed = all(n_raw_cl <= maxClosed) && all(n_cln_cl <= maxClosed);
                summary.valid_open   = all(n_raw_op <= maxOpen)   && all(n_cln_op <= maxOpen);
                summary.maxBlinkRate = maxBlinkRate;
                summary.amplified_closed = BlinkAnalysis.blinkAmplification(n_raw_cl, n_cln_cl);
                summary.amplified_open   = BlinkAnalysis.blinkAmplification(n_raw_op, n_cln_op);

                obj.results(c).perChannel = perCh;
                obj.results(c).summary    = summary;
                obj.results(c).method     = method;
            end

            fprintf('  Done.\n');
        end

        % ── getSummary ───────────────────────────────────────────────
        function vals = getSummary(obj, fieldName)
            % GETSUMMARY - Pull one summary scalar across all combos.
            %   Returns [1 x nCombos] vector.
            obj.checkComputed();
            vals = zeros(1, obj.nCombos);
            for c = 1:obj.nCombos
                vals(c) = obj.results(c).summary.(fieldName);
            end
        end

        % ── paramTable ───────────────────────────────────────────────
        function T = paramTable(obj)
            % PARAMTABLE - One row per combo: parameters + key blink scalars.
            obj.checkComputed();
            rows = cell(obj.nCombos, 1);
            for c = 1:obj.nCombos
                params = obj.subject_results(c).parameters;
                fields = fieldnames(params);
                row    = table();
                row.comboIdx = c;
                for fi = 1:numel(fields)
                    val = params.(fields{fi});
                    if isscalar(val) && isnumeric(val)
                        row.(fields{fi}) = val;
                    end
                end
                s = obj.results(c).summary;
                row.reductionRatio_closed  = s.blinkReductionRatio_closed;
                row.reductionPct_closed    = s.blinkReductionPct_closed;
                row.persistenceRate_closed = s.persistenceRate_closed;
                row.removalRate_closed     = s.removalRate_closed;
                row.artefactRate_closed    = s.artefactRate_closed;
                row.reductionRatio_open    = s.blinkReductionRatio_open;
                row.persistenceRate_open   = s.persistenceRate_open;
                row.removalRate_open       = s.removalRate_open;
                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function checkComputed(obj)
            if isempty(obj.results) || ~isfield(obj.results(1), 'summary')
                error('BlinkAnalysis:notComputed', 'Call compute() before accessing results.');
            end
        end

    end

    % ================================================================
    methods (Static, Access = public)

        % ── detectBlinks ─────────────────────────────────────────────
        function [n_raw, n_clean, idx_raw, idx_clean] = ...
                detectBlinks(XRaw, XClean, thresh, minDistRaw, minDistClean, fs, b, a, method)
            % DETECTBLINKS - Per-channel peak detection with asymmetric minPeakDistance.

            nCh  = size(XRaw, 1);
            mpdR = round(minDistRaw   * fs);
            mpdC = round(minDistClean * fs);

            idx_raw   = cell(nCh, 1);
            idx_clean = cell(nCh, 1);

            for ch = 1:nCh
                if strcmp(method, 'mad')
                    xr = filtfilt(b, a, double(XRaw(ch,:)));
                    xc = filtfilt(b, a, double(XClean(ch,:)));
                else
                    xr = double(XRaw(ch,:));
                    xc = double(XClean(ch,:));
                end
                [~, idx_raw{ch}]   = findpeaks(abs(xr), ...
                    'MinPeakHeight', thresh(ch), 'MinPeakDistance', mpdR);
                [~, idx_clean{ch}] = findpeaks(abs(xc), ...
                    'MinPeakHeight', thresh(ch), 'MinPeakDistance', mpdC);
            end

            n_raw   = cellfun(@numel, idx_raw);
            n_clean = cellfun(@numel, idx_clean);
        end

        % ── matchBlinks ──────────────────────────────────────────────
        function [nMatched, idxMatchedRaw, idxMatchedClean] = ...
                matchBlinks(idx_raw_cells, idx_clean_cells, tolSamples, offset)
            % MATCHBLINKS - Greedy 1-D bipartite matching for blink persistence.
            %
            %   Each raw blink and each clean blink is used at most once.
            %   Guarantees: nMatched <= min(n_raw, n_clean).
            %
            %   Algorithm:
            %     For each clean blink (sorted), find the nearest unclaimed
            %     raw blink within ±tolSamples.  Claim it.  Repeat.
            %
            %   offset : added to clean indices before matching (0 for our
            %            pipeline; non-zero if clean was trimmed externally).
            %            Matched clean indices stored WITHOUT offset (original
            %            space) so callers can overlay them on the clean signal.

            nCh             = numel(idx_raw_cells);
            nMatched        = zeros(nCh, 1);
            idxMatchedRaw   = cell(nCh, 1);
            idxMatchedClean = cell(nCh, 1);

            for ch = 1:nCh
                rawIdx   = sort(double(idx_raw_cells{ch}(:)));
                cleanIdx = sort(double(idx_clean_cells{ch}(:))) + offset;

                if isempty(rawIdx) || isempty(cleanIdx)
                    idxMatchedRaw{ch}   = zeros(0,1);
                    idxMatchedClean{ch} = zeros(0,1);
                    continue;
                end

                claimed      = false(numel(rawIdx), 1);
                matchedRaw   = zeros(0,1);
                matchedClean = zeros(0,1);

                for ci = 1:numel(cleanIdx)
                    cIdx = cleanIdx(ci);

                    % Distance to each raw blink; Inf for already-claimed ones
                    dists         = abs(rawIdx - cIdx);
                    dists(claimed) = Inf;

                    % Best candidate must be within tolerance
                    [bestDist, bestRaw] = min(dists);
                    if bestDist > tolSamples
                        continue;   % no valid match for this clean blink
                    end

                    % Claim
                    claimed(bestRaw)       = true;
                    matchedRaw(end+1,1)    = rawIdx(bestRaw);
                    matchedClean(end+1,1)  = cleanIdx(ci) - offset;  % store in original space
                end

                nMatched(ch)          = numel(matchedRaw);
                idxMatchedRaw{ch}     = matchedRaw;
                idxMatchedClean{ch}   = matchedClean;
            end
        end

        % ── Rate helper functions ─────────────────────────────────────
        function r = reductionRatio(n_raw, n_clean)
            total_raw = sum(n_raw);
            if total_raw == 0, r = NaN;
            else, r = max(0, 1 - sum(n_clean) / total_raw); end
        end

        function r = persistenceRate(n_raw, n_matched)
            % Fraction of raw blinks NOT removed (0=perfect, 1=ASR did nothing)
            total_raw = sum(n_raw);
            if total_raw == 0, r = NaN;
            else, r = sum(n_matched) / total_raw; end
        end

        function r = removalRate(n_raw, n_matched)
            % Fraction of raw blinks successfully removed
            total_raw = sum(n_raw);
            if total_raw == 0, r = NaN;
            else, r = 1 - sum(n_matched) / total_raw; end
        end

        function r = artefactRate(n_clean, n_matched)
            % Fraction of clean blinks that are NEW (not matching any raw blink)
            total_clean = sum(n_clean);
            if total_clean == 0, r = NaN;
            else, r = max(0, 1 - sum(n_matched) / total_clean); end
        end

        function flag = blinkAmplification(n_raw, n_clean)
            flag = sum(n_clean) > sum(n_raw);
        end

    end % static methods

end
