classdef ICAAnalysis < handle
    % ICAAnalysis - ICA-based evaluation of ASR cleaning quality.
    %
    %   Implements the analysis framework of Chang et al. (2020), IEEE TBME:
    %   "Evaluation of Artifact Subspace Reconstruction for Automatic Artifact
    %   Components Removal in Multi-Channel EEG Recordings."
    %
    %   The core idea: run ICA on both the raw and ASR-cleaned signals, then
    %   quantify what ASR did in terms of independent components rather than
    %   raw signal amplitude. This separates artifact removal (good) from brain
    %   signal removal (bad) in a way that sample-level metrics like RRMSE
    %   cannot.
    %
    %   ── THREE MEASURES (directly from Chang et al.) ──────────────────────
    %
    %   1. IC PRESERVATION (Section III-B, Fig. 3/4)
    %      Matching ICs between raw and cleaned ICA decompositions via the
    %      Hungarian algorithm. ICs with correlation > preservationThresh
    %      (default 0.8) are "preserved". Reported per class (Brain, Eye,
    %      Muscle, Other).
    %
    %   2. SOURCE POWER REDUCTION (Section III-C, Fig. 5)
    %      Apply the RAW unmixing matrix W to the CLEANED data (Eq. 2 in
    %      paper): Y(k) = W_raw * X_clean. Power retained per IC =
    %      Var(Y_clean) / Var(Y_raw). This is the key metric because it uses
    %      a fixed spatial filter (from raw), so changes reflect what ASR
    %      actually removed from each source, not a different decomposition.
    %
    %   3. IC CLASSIFICATION via ICLabel (Section II-C3)
    %      Classifies each IC into: Brain, Eye, Muscle, Heart, Line Noise,
    %      Channel Noise, Other. Requires EEGLAB's iclabel plugin.
    %      Gracefully skips if iclabel is not available (labels all as 'Other').
    %
    %   4. DIPOLE FITTING (Section II-C4) — OPTIONAL
    %      Requires dipfit plugin and 3D chanlocs. Enabled by calling
    %      setDipfitParams(hdmFile, mriFile, chanlocs) before compute().
    %      If not configured, dipolarity is skipped silently.
    %
    %   ── EEGLAB DEPENDENCY ────────────────────────────────────────────────
    %   Requires EEGLAB in the MATLAB path. Functions used:
    %     runica()      — extended Infomax ICA (always needed)
    %     iclabel()     — IC classification (graceful skip if absent)
    %     pop_dipfit*   — dipole fitting (optional, skip if absent)
    %   Check with: which('runica') ~= 0
    %
    %   ── CONSTRUCTOR ──────────────────────────────────────────────────────
    %   Identical signature to other analyser classes:
    %       ica = ICAAnalysis(raw, subject_results, fs, subjectID, algName);
    %       ica = ICAAnalysis(raw, subject_results, fs, subjectID, algName, chanlocs);
    %
    %   ── COMPUTE ──────────────────────────────────────────────────────────
    %       ica.compute();                         % all combos, default options
    %       ica.compute('segment', 'closed+open',  % which data to decompose
    %                   'maxSamples', 60000,        % cap for speed (default 60000)
    %                   'preservationThresh', 0.8); % IC match threshold
    %
    %   ── RESULTS ──────────────────────────────────────────────────────────
    %       ica.results(c).raw              — ICA of raw signal
    %       ica.results(c).cleaned          — ICA of cleaned signal
    %       ica.results(c).comparison       — matching + power metrics
    %       ica.results(c).summary          — scalar summary per combo
    %
    %       ica.getSummary('powerRetained_brain')   % [1 x nCombos]
    %       ica.getSummary('preservationRate_eye')
    %       ica.paramTable()
    %
    %   ── INTEGRATION WITH ExperimentAnalysis ──────────────────────────────
    %       ana.computeICA();               % adds ana.ica
    %       ana.ica.results(k).comparison.powerRetained_byClass
    %
    %   ── COMPUTATIONAL COST WARNING ───────────────────────────────────────
    %   runica on 32-ch x 60000 samples takes ~2-10 minutes per run.
    %   With 10 combos x 2 ICA decompositions per combo = up to 200 ICA runs.
    %   STRONGLY RECOMMENDED: set maxSamples (default 60000) and run on HPC.
    %   Alternatively, run ICA on raw once and reuse W for cleaned (via
    %   'reuseRawW', true option) — this skips the cleaned ICA decomposition.
    %
    %   Reference:
    %     Chang et al. (2020). Evaluation of Artifact Subspace Reconstruction
    %     for Automatic Artifact Components Removal in Multi-Channel EEG
    %     Recordings. IEEE Trans. Biomed. Eng., 67(4), 1114-1121.

    % ================================================================
    properties (Access = public)

        % ── Inputs ───────────────────────────────────────────────────
        raw              % struct: .calibration .closed .open  [C x T each]
        subject_results  % 1xN struct (same layout as other analysers)
        fs               % sampling rate Hz
        subjectID
        algorithmName
        chanlocs         % EEGLAB chanlocs struct (optional — improves ICLabel)

        % ── Results ──────────────────────────────────────────────────
        results          % 1 x nCombos struct array — see header for fields
        nCombos

        % ── Parameters (set by compute()) ────────────────────────────
        segment              = 'closed+open'   % which segments to decompose
        maxSamples           = 60000           % max samples per ICA (cap for speed)
        preservationThresh   = 0.8             % IC correlation threshold for "preserved"
        reuseRawW            = false           % skip cleaned ICA — apply raw W directly
        icLabelThresh        = 0.7             % min probability to assign a class label

        % ── Dipfit (optional) ────────────────────────────────────────
        dipfitEnabled  = false
        dipfitHDM      = ''    % path to head model file (.mat)
        dipfitMRI      = ''    % path to MRI file (e.g. standard_BEM)
        dipfitChanlocs = []    % 3D chanlocs struct for dipfit

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = ICAAnalysis(raw, subject_results, fs, subjectID, algorithmName, chanlocs)
            % ICAANALYSIS
            %   raw            : struct with .calibration/.closed/.open
            %   subject_results: 1xN struct array with .cleanClosed/.cleanOpen etc.
            %   fs             : sampling rate Hz
            %   subjectID      : numeric subject identifier
            %   algorithmName  : string (e.g. 'ema-asr')
            %   chanlocs       : (optional) EEGLAB chanlocs struct

            obj.raw             = raw;
            obj.subject_results = subject_results;
            obj.fs              = fs;
            obj.subjectID       = subjectID;
            obj.algorithmName   = algorithmName;
            obj.nCombos         = numel(subject_results);

            if nargin >= 6 && ~isempty(chanlocs)
                obj.chanlocs = chanlocs;
            end

            obj.checkEEGLAB();
        end

        % ── setDipfitParams ───────────────────────────────────────────
        function setDipfitParams(obj, hdmFile, mriFile, chanlocs3D)
            % SETDIPFITPARAMS - Enable dipole fitting.
            %   hdmFile    : path to head model .mat (e.g. dipfit standard_BEM)
            %   mriFile    : path to MRI .mat
            %   chanlocs3D : EEGLAB chanlocs struct with 3D coordinates
            %
            %   Example (using EEGLAB defaults):
            %     eeglab_path = fileparts(which('eeglab'));
            %     hdm = fullfile(eeglab_path,'plugins','dipfit','standard_BEM','standard_vol.mat');
            %     mri = fullfile(eeglab_path,'plugins','dipfit','standard_BEM','standard_mri.mat');
            %     ica.setDipfitParams(hdm, mri, EEG.chanlocs);

            if ~isfile(hdmFile)
                error('ICAAnalysis:dipfit', 'Head model file not found: %s', hdmFile);
            end
            obj.dipfitHDM      = hdmFile;
            obj.dipfitMRI      = mriFile;
            obj.dipfitChanlocs = chanlocs3D;
            obj.dipfitEnabled  = true;
            fprintf('[ICAAnalysis] Dipfit enabled.\n');
        end

        % ── compute ──────────────────────────────────────────────────
        function compute(obj, varargin)
            % COMPUTE - Run ICA analysis across all parameter combinations.
            %
            %   Name-value options:
            %     'segment'            'closed' | 'open' | 'closed+open' (default)
            %                          | 'all'  — which segments to decompose
            %     'maxSamples'         scalar (default 60000) — cap on samples
            %                          fed to runica. Data is randomly subsampled
            %                          if longer. Set Inf to use all samples.
            %     'preservationThresh' scalar (default 0.8) — correlation cutoff
            %                          for classifying an IC as "preserved"
            %     'reuseRawW'          logical (default false) — if true, skip the
            %                          separate ICA on cleaned data and apply the
            %                          raw W directly. Cuts compute time in half.
            %     'icLabelThresh'      scalar (default 0.7) — minimum ICLabel
            %                          probability to assign a class (else 'Other')

            p = inputParser;
            p.KeepUnmatched = false;
            addParameter(p, 'segment',            obj.segment,            @ischar);
            addParameter(p, 'maxSamples',         obj.maxSamples,         @isnumeric);
            addParameter(p, 'preservationThresh', obj.preservationThresh, @isnumeric);
            addParameter(p, 'reuseRawW',          obj.reuseRawW,          @(x) islogical(x)||isnumeric(x));
            addParameter(p, 'icLabelThresh',      obj.icLabelThresh,      @isnumeric);
            parse(p, varargin{:});

            obj.segment            = p.Results.segment;
            obj.maxSamples         = p.Results.maxSamples;
            obj.preservationThresh = p.Results.preservationThresh;
            obj.reuseRawW          = logical(p.Results.reuseRawW);
            obj.icLabelThresh      = p.Results.icLabelThresh;

            fprintf('\n=== ICAAnalysis: S%d [%s] ===\n', obj.subjectID, obj.algorithmName);
            fprintf('  Combos    : %d\n',    obj.nCombos);
            fprintf('  Segment   : %s\n',    obj.segment);
            fprintf('  MaxSamples: %d\n',    obj.maxSamples);
            fprintf('  ReuseRawW : %d\n\n',  obj.reuseRawW);

            % ── Assemble raw data once — shared across all combos ─────
            X_raw = obj.assembleSegment(obj.raw.closed, obj.raw.open, ...
                                        obj.raw.calibration, obj.segment);
            X_raw = obj.subsampleData(X_raw);

            fprintf('  Running ICA on raw signal [%d ch x %d samples]...\n', ...
                size(X_raw, 1), size(X_raw, 2));
            [W_raw, A_raw, act_raw] = obj.runICA(X_raw);
            icLabels_raw = obj.classifyICs(W_raw, A_raw, act_raw, X_raw);

            fprintf('  Raw ICA done. Running per-combo analysis...\n\n');

            % ── Per-combo loop ────────────────────────────────────────
            for c = 1:obj.nCombos
                sr = obj.subject_results(c);
                fprintf('  [%d/%d] combo %d\n', c, obj.nCombos, c);

                % Assemble cleaned data for this combo
                X_clean = obj.assembleSegment( ...
                    sr.cleanClosed, sr.cleanOpen, sr.cleanCalibration, obj.segment);
                X_clean = obj.subsampleData(X_clean);

                % ICA on cleaned (or reuse raw W)
                if obj.reuseRawW
                    W_clean  = W_raw;
                    A_clean  = A_raw;
                    act_clean = W_clean * X_clean;
                    icLabels_clean = icLabels_raw;   % same topology assumed
                    fprintf('    (reusing raw W for cleaned)\n');
                else
                    fprintf('    Running ICA on cleaned [%d ch x %d samples]...\n', ...
                        size(X_clean, 1), size(X_clean, 2));
                    [W_clean, A_clean, act_clean] = obj.runICA(X_clean);
                    icLabels_clean = obj.classifyICs(W_clean, A_clean, act_clean, X_clean);
                end

                % Apply raw W to cleaned data — Eq. 2 in Chang et al.
                act_raw_on_clean = W_raw * X_clean;

                % IC preservation via Hungarian matching
                [matchedPairs, icCorr] = obj.matchICs(A_raw, A_clean);
                corrDiag   = diag(icCorr(matchedPairs(:,1), matchedPairs(:,2)));
                preserved  = corrDiag >= obj.preservationThresh;

                % Per-IC source power retention
                [powerRetained, varRaw, varClean] = ...
                    obj.computePowerRetained(act_raw, act_raw_on_clean);

                % Class-wise preservation and power
                classNames = {'Brain','Eye','Muscle','Heart','Line Noise','Channel Noise','Other'};
                compOut    = obj.buildComparison( ...
                    matchedPairs, icCorr, corrDiag, preserved, ...
                    powerRetained, varRaw, varClean, ...
                    icLabels_raw, icLabels_clean, classNames);

                % Dipole fitting (optional)
                dipolarity_raw   = [];
                dipolarity_clean = [];
                if obj.dipfitEnabled
                    dipolarity_raw   = obj.runDipfit(W_raw,   A_raw,   X_raw);
                    dipolarity_clean = obj.runDipfit(W_clean, A_clean, X_clean);
                end

                % Pack results
                obj.results(c).raw = struct( ...
                    'W',            W_raw,           ...
                    'A',            A_raw,           ...
                    'icaact',       act_raw,         ...
                    'iclabels',     icLabels_raw.probs,  ...
                    'icClass',      {icLabels_raw.classes}, ...
                    'varExplained', obj.varianceExplained(A_raw, X_raw), ...
                    'dipolarity',   dipolarity_raw);

                obj.results(c).cleaned = struct( ...
                    'W',            W_clean,         ...
                    'A',            A_clean,         ...
                    'icaact',       act_clean,       ...
                    'iclabels',     icLabels_clean.probs, ...
                    'icClass',      {icLabels_clean.classes}, ...
                    'varExplained', obj.varianceExplained(A_clean, X_clean), ...
                    'dipolarity',   dipolarity_clean);

                obj.results(c).comparison = compOut;
                obj.results(c).summary    = obj.buildSummary(compOut, ...
                    dipolarity_raw, dipolarity_clean, icLabels_raw);
                obj.results(c).parameters = sr.parameters;

                fprintf('    Preserved: %d/%d ICs (corr>%.1f)  Brain power retained: %.0f%%\n', ...
                    sum(preserved), numel(preserved), obj.preservationThresh, ...
                    obj.results(c).summary.powerRetained_brain * 100);
            end

            fprintf('\n=== ICAAnalysis complete ===\n');
        end

        % ── getSummary ────────────────────────────────────────────────
        function vals = getSummary(obj, fieldName)
            % GETSUMMARY - Extract one scalar summary field across all combos.
            %   fieldName : e.g. 'powerRetained_brain', 'preservationRate',
            %               'preservationRate_eye', 'nICsPreserved',
            %               'nDipolarBrain_raw', 'nDipolarBrain_clean'
            if isempty(obj.results)
                error('ICAAnalysis:notComputed', 'Call compute() first.');
            end
            vals = zeros(1, obj.nCombos);
            for c = 1:obj.nCombos
                s = obj.results(c).summary;
                if isfield(s, fieldName)
                    v = s.(fieldName);
                    if isnumeric(v) && isscalar(v)
                        vals(c) = v;
                    else
                        vals(c) = NaN;
                    end
                else
                    vals(c) = NaN;
                end
            end
        end

        % ── paramTable ────────────────────────────────────────────────
        function T = paramTable(obj, segment)
            % PARAMTABLE - One row per combo: parameters + key ICA metrics.
            if nargin < 2, segment = obj.segment; end
            if isempty(obj.results)
                error('ICAAnalysis:notComputed', 'Call compute() first.');
            end

            rows = cell(obj.nCombos, 1);
            for c = 1:obj.nCombos
                params = obj.results(c).parameters;
                s      = obj.results(c).summary;
                row    = table();
                row.comboIdx = c;

                % Parameter columns
                f = fieldnames(params);
                for i = 1:numel(f)
                    v = params.(f{i});
                    if isscalar(v) && isnumeric(v)
                        row.(f{i}) = v;
                    end
                end

                % Key metric columns
                metrics_to_add = { ...
                    'nICsPreserved', 'preservationRate', ...
                    'preservationRate_brain', 'preservationRate_eye', ...
                    'preservationRate_muscle', ...
                    'powerRetained_brain', 'powerRetained_eye', ...
                    'powerRetained_muscle', ...
                    'nDipolarBrain_raw', 'nDipolarBrain_clean'};

                for i = 1:numel(metrics_to_add)
                    fn = metrics_to_add{i};
                    if isfield(s, fn) && ~isempty(s.(fn)) && isnumeric(s.(fn))
                        row.(fn) = s.(fn);
                    end
                end

                rows{c} = row;
            end
            T = vertcat(rows{:});
        end

        % ── plotPreservation ─────────────────────────────────────────
        function plotPreservation(obj)
            % PLOTPRESERVATION - Bar chart of IC preservation rate by class.
            if isempty(obj.results)
                error('ICAAnalysis:notComputed', 'Call compute() first.');
            end

            classes = {'Brain', 'Eye', 'Muscle', 'Other'};
            fields  = {'preservationRate_brain', 'preservationRate_eye', ...
                       'preservationRate_muscle', 'preservationRate_other'};
            nC      = obj.nCombos;
            data    = zeros(4, nC);

            for c = 1:nC
                s = obj.results(c).summary;
                for i = 1:4
                    if isfield(s, fields{i})
                        data(i,c) = s.(fields{i});
                    end
                end
            end

            figure('Name', sprintf('IC Preservation — S%d %s', ...
                obj.subjectID, obj.algorithmName));
            bar(data' * 100);
            legend(classes, 'Location', 'best');
            xlabel('Parameter combo'); ylabel('Preserved ICs (%)');
            title(sprintf('IC Preservation Rate — S%d [%s]', ...
                obj.subjectID, obj.algorithmName));
            ylim([0 105]);
            grid on;
        end

        % ── plotPowerRetained ─────────────────────────────────────────
        function plotPowerRetained(obj)
            % PLOTPOWERRETAINED - Power retention by IC class across combos.
            if isempty(obj.results)
                error('ICAAnalysis:notComputed', 'Call compute() first.');
            end

            classes = {'Brain', 'Eye', 'Muscle'};
            fields  = {'powerRetained_brain', 'powerRetained_eye', 'powerRetained_muscle'};
            nC      = obj.nCombos;
            data    = zeros(3, nC);

            for c = 1:nC
                s = obj.results(c).summary;
                for i = 1:3
                    if isfield(s, fields{i})
                        data(i,c) = s.(fields{i});
                    end
                end
            end

            figure('Name', sprintf('Power Retained — S%d %s', ...
                obj.subjectID, obj.algorithmName));
            plot(data' * 100, 'o-', 'LineWidth', 1.5);
            legend(classes, 'Location', 'best');
            xlabel('Parameter combo'); ylabel('Power retained (%)');
            title(sprintf('Source Power Retained — S%d [%s]', ...
                obj.subjectID, obj.algorithmName));
            ylim([0 105]); grid on;
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        % ── assembleSegment ───────────────────────────────────────────
        function X = assembleSegment(obj, closed, open, calibration, seg)
            switch lower(seg)
                case 'closed',      X = closed;
                case 'open',        X = open;
                case 'closed+open', X = [closed, open];
                case 'all',         X = [calibration, closed, open];
                otherwise
                    error('ICAAnalysis:badSegment', ...
                        'Unknown segment "%s". Use: closed|open|closed+open|all', seg);
            end
            X = double(X);
            X(~isfinite(X)) = 0;
        end

        % ── subsampleData ─────────────────────────────────────────────
        function X = subsampleData(obj, X)
            % SUBSAMPLEDATA - Random column subsample to maxSamples.
            %   Preserves temporal order (sorted random indices).
            [~, T] = size(X);
            if T > obj.maxSamples && isfinite(obj.maxSamples)
                idx = sort(randperm(T, round(obj.maxSamples)));
                X   = X(:, idx);
            end
        end

        % ── runICA ────────────────────────────────────────────────────
        function [W, A, icaact] = runICA(~, X)
            % RUNICA - Extended Infomax ICA via EEGLAB's runica.
            %   Returns:
            %     W      : unmixing matrix [C x C]
            %     A      : mixing matrix   [C x C] = pinv(W)
            %     icaact : IC activations  [C x T] = W * X

            % Remove mean per channel (ICA requires zero-mean)
            X = X - mean(X, 2);

            % runica returns weights + sphere (W = weights * sphere)
            [weights, sphere] = runica(X, 'extended', 1, 'verbose', 'off');
            W      = weights * sphere;
            A      = pinv(W);
            icaact = W * X;
        end

        % ── classifyICs ───────────────────────────────────────────────
        function labels = classifyICs(obj, W, A, icaact, X)
            % CLASSIFYICS - Classify ICs using ICLabel.
            %   Returns struct with fields:
            %     .probs   [C x 7] probability per class
            %     .classes {C x 1} cell array of string class labels

            classNames = {'Brain','Eye','Muscle','Heart','Line Noise','Channel Noise','Other'};
            C = size(W, 1);
            T = size(X,  2);

            % Build minimal EEGLAB EEG struct
            EEG              = struct();
            EEG.data         = X;
            EEG.icaweights   = W;
            EEG.icasphere    = eye(C);
            EEG.icawinv      = A;
            EEG.icaact       = icaact;
            EEG.srate        = obj.fs;
            EEG.nbchan       = C;
            EEG.trials       = 1;
            EEG.pnts         = T;
            EEG.xmin         = 0;
            EEG.times        = (0:T-1) / obj.fs;
            EEG.chanlocs     = obj.chanlocs;     % [] is fine — ICLabel degrades gracefully
            EEG.chaninfo     = struct();
            EEG.ref          = 'averef';
            EEG.event        = [];
            EEG.epoch        = [];

            % Run ICLabel
            try
                EEG     = iclabel(EEG);
                probs   = EEG.etc.ic_classification.ICLabel.classifications;   % [C x 7]
            catch e
                warning('ICAAnalysis:iclabelFailed', ...
                    'iclabel() failed: %s\nLabelling all ICs as Other.', e.message);
                probs = zeros(C, 7);
                probs(:, 7) = 1;   % all → Other
            end

            % Assign class label by argmax, subject to icLabelThresh
            [maxProb, maxIdx] = max(probs, [], 2);
            classes = repmat({'Other'}, C, 1);
            for i = 1:C
                if maxProb(i) >= obj.icLabelThresh
                    classes{i} = classNames{maxIdx(i)};
                end
            end

            labels.probs   = probs;
            labels.classes = classes;
        end

        % ── matchICs ─────────────────────────────────────────────────
        function [matchedPairs, corrMat] = matchICs(~, A_raw, A_clean)
            % MATCHICS - Hungarian-algorithm matching of ICs by mixing column correlation.
            %   Matches columns of A_raw to columns of A_clean to maximise
            %   the sum of absolute correlations (as in Chang et al., Eq. 3).
            %
            %   Returns:
            %     matchedPairs : [C x 2] each row is [raw_ic_idx, clean_ic_idx]
            %     corrMat      : [C x C] full correlation matrix (raw cols vs clean cols)

            C = size(A_raw, 2);

            % Pearson correlation between all column pairs
            corrMat = zeros(C);
            for i = 1:C
                for j = 1:C
                    a = A_raw(:,i)   - mean(A_raw(:,i));
                    b = A_clean(:,j) - mean(A_clean(:,j));
                    na = norm(a); nb = norm(b);
                    if na > 0 && nb > 0
                        corrMat(i,j) = abs(dot(a,b) / (na*nb));
                    end
                end
            end

            % Hungarian assignment to maximise total correlation
            % MATLAB's matchpairs() minimises cost, so use 1 - corrMat as cost.
            costMat = 1 - corrMat;
            try
                % matchpairs available in R2020b+
                matchedPairs = matchpairs(costMat, 1.01);   % padCost > max(costMat)
            catch
                % Fallback: greedy matching (good enough for small C)
                matchedPairs = ICAAnalysis.greedyMatch(corrMat);
            end
        end

        % ── computePowerRetained ──────────────────────────────────────
        function [powerRetained, varRaw, varClean] = computePowerRetained(~, act_raw, act_raw_on_clean)
            % COMPUTEPOWERRETAINED - Per-IC source power retention (Chang et al. Eq. 2-3).
            %   act_raw          : [C x T] IC activations of raw data  (W_raw * X_raw)
            %   act_raw_on_clean : [C x T] IC activations of cleaned   (W_raw * X_clean)
            %
            %   powerRetained(i) = Var(act_raw_on_clean(i,:)) / Var(act_raw(i,:))

            varRaw   = var(act_raw,           0, 2);   % [C x 1]
            varClean = var(act_raw_on_clean,  0, 2);   % [C x 1]

            % Guard against zero-variance ICs
            safeVar          = max(varRaw, eps);
            powerRetained    = varClean ./ safeVar;
            powerRetained    = min(powerRetained, 1);  % cap at 1 (ASR can't add power)
        end

        % ── buildComparison ───────────────────────────────────────────
        function comp = buildComparison(obj, matchedPairs, icCorr, corrDiag, ...
                                        preserved, powerRetained, varRaw, varClean, ...
                                        icLabels_raw, icLabels_clean, classNames)

            C = size(icCorr, 1);

            % Preservation by class
            rawClasses = icLabels_raw.classes;   % {C x 1}
            classMap   = struct('Brain', 1, 'Eye', 2, 'Muscle', 3, 'Other', 7);
            mainClasses = {'Brain', 'Eye', 'Muscle', 'Other'};

            preservedByClass  = struct();
            powerByClass      = struct();
            classCountsRaw    = struct();
            classCountsClean  = struct();

            for ci = 1:numel(mainClasses)
                cls = mainClasses{ci};
                fn  = lower(cls);

                % Which raw ICs are this class?
                isClass = strcmp(rawClasses, cls);
                % For 'Other': aggregate all non-Brain/Eye/Muscle/Heart/Line/Channel
                if strcmp(cls, 'Other')
                    isClass = ~strcmp(rawClasses, 'Brain') & ...
                              ~strcmp(rawClasses, 'Eye')   & ...
                              ~strcmp(rawClasses, 'Muscle');
                end

                classCountsRaw.(fn) = sum(isClass);

                if sum(isClass) > 0
                    % Preservation: matched raw ICs of this class that are preserved
                    matchedRawIdx = matchedPairs(:,1);
                    classMatchMask = isClass(matchedRawIdx);
                    preservedByClass.(fn) = sum(preserved(classMatchMask)) / ...
                                            max(sum(classMatchMask), 1);
                    powerByClass.(fn) = mean(powerRetained(isClass));
                else
                    preservedByClass.(fn) = NaN;
                    powerByClass.(fn)     = NaN;
                end

                % Count in clean decomposition
                isClassClean = strcmp(icLabels_clean.classes, cls);
                if strcmp(cls, 'Other')
                    isClassClean = ~strcmp(icLabels_clean.classes, 'Brain') & ...
                                   ~strcmp(icLabels_clean.classes, 'Eye')   & ...
                                   ~strcmp(icLabels_clean.classes, 'Muscle');
                end
                classCountsClean.(fn) = sum(isClassClean);
            end

            comp.icCorr              = icCorr;
            comp.matchedPairs        = matchedPairs;
            comp.corrPerMatchedPair  = corrDiag;
            comp.preservedMask       = preserved;
            comp.preservedIdx        = find(preserved);
            comp.preservation        = mean(preserved);
            comp.preservedByClass    = preservedByClass;
            comp.powerRetained       = powerRetained;
            comp.powerByClass        = powerByClass;
            comp.varRaw              = varRaw;
            comp.varClean            = varClean;
            comp.classCountsRaw      = classCountsRaw;
            comp.classCountsClean    = classCountsClean;
        end

        % ── buildSummary ──────────────────────────────────────────────
        function s = buildSummary(obj, comp, dipolarity_raw, dipolarity_clean, icLabels_raw)
            s.nICsPreserved      = sum(comp.preservedMask);
            s.preservationRate   = comp.preservation;

            % Per-class preservation rates (flat fields for easy getSummary access)
            classes = {'brain', 'eye', 'muscle', 'other'};
            for ci = 1:numel(classes)
                cls = classes{ci};
                fn  = ['preservationRate_', cls];
                if isfield(comp.preservedByClass, cls)
                    s.(fn) = comp.preservedByClass.(cls);
                else
                    s.(fn) = NaN;
                end
            end

            % Per-class power retained
            for ci = 1:numel(classes)
                cls = classes{ci};
                fn  = ['powerRetained_', cls];
                if isfield(comp.powerByClass, cls)
                    s.(fn) = comp.powerByClass.(cls);
                else
                    s.(fn) = NaN;
                end
            end

            % Dipole fitting results (optional)
            s.nDipolarBrain_raw   = NaN;
            s.nDipolarBrain_clean = NaN;
            if ~isempty(dipolarity_raw)
                isBrain = strcmp(icLabels_raw.classes, 'Brain');
                s.nDipolarBrain_raw   = sum(dipolarity_raw(isBrain));
            end
            if ~isempty(dipolarity_clean)
                s.nDipolarBrain_clean = sum(dipolarity_clean);
            end
        end

        % ── varianceExplained ─────────────────────────────────────────
        function ve = varianceExplained(~, A, X)
            % VARIANCEEXPLAINED - % variance explained by each IC.
            %   Uses the standard EEG formula: project IC back to scalp,
            %   compute its variance as fraction of total signal variance.
            C      = size(A, 2);
            icaact = pinv(A) * X;
            totalVar = sum(var(X, 0, 2));
            ve     = zeros(C, 1);
            for i = 1:C
                reconst  = A(:,i) * icaact(i,:);
                ve(i)    = sum(var(reconst, 0, 2)) / max(totalVar, eps);
            end
        end

        % ── runDipfit ─────────────────────────────────────────────────
        function dipolarity = runDipfit(obj, W, A, X)
            % RUNDIPFIT - Dipole fitting via EEGLAB dipfit.
            %   Returns [C x 1] logical: true if residual variance < 5%.
            %   Requires obj.dipfitEnabled = true and valid HDM/MRI paths.

            C = size(W, 1);
            T = size(X, 2);
            dipolarity = false(C, 1);

            if ~obj.dipfitEnabled
                return;
            end

            try
                EEG              = struct();
                EEG.data         = X;
                EEG.icaweights   = W;
                EEG.icasphere    = eye(C);
                EEG.icawinv      = A;
                EEG.icaact       = W * X;
                EEG.srate        = obj.fs;
                EEG.nbchan       = C;
                EEG.trials       = 1;
                EEG.pnts         = T;
                EEG.xmin         = 0;
                EEG.times        = (0:T-1) / obj.fs;
                EEG.chanlocs     = obj.dipfitChanlocs;
                EEG.chaninfo     = struct();
                EEG.ref          = 'averef';
                EEG.event        = [];
                EEG.epoch        = [];

                % Set dipfit model
                EEG = pop_dipfit_settings(EEG, ...
                    'hdmfile',  obj.dipfitHDM, ...
                    'mrifile',  obj.dipfitMRI, ...
                    'chanfile', '', ...
                    'coordformat', 'MNI');

                % Fit dipoles
                EEG = pop_multifit(EEG, 1:C, 'threshold', 100, 'dipoles', 1);

                % Extract residual variance
                for i = 1:C
                    if ~isempty(EEG.dipfit.model(i).rv)
                        dipolarity(i) = EEG.dipfit.model(i).rv < 0.05;
                    end
                end

            catch e
                warning('ICAAnalysis:dipfitFailed', ...
                    'Dipfit failed: %s\nSkipping dipolarity.', e.message);
                dipolarity = false(C, 1);
            end
        end

        % ── checkEEGLAB ───────────────────────────────────────────────
        function checkEEGLAB(~)
            if isempty(which('runica'))
                error('ICAAnalysis:noEEGLAB', ...
                    ['EEGLAB not found in MATLAB path.\n' ...
                     'Add EEGLAB to the path before using ICAAnalysis:\n' ...
                     '  addpath(genpath(''/path/to/eeglab''));\n' ...
                     '  eeglab nogui;']);
            end
            if isempty(which('iclabel'))
                warning('ICAAnalysis:noICLabel', ...
                    ['iclabel plugin not found. IC classification will label\n' ...
                     'all ICs as Other. Install via EEGLAB plugin manager:\n' ...
                     '  File > Manage EEGLAB extensions > ICLabel']);
            end
        end

    end % private methods

    % ================================================================
    methods (Static, Access = private)

        function pairs = greedyMatch(corrMat)
            % GREEDYMATCH - Greedy fallback when matchpairs() is unavailable.
            %   Repeatedly picks the highest unclaimed correlation.
            C     = size(corrMat, 1);
            pairs = zeros(C, 2);
            used_raw   = false(C, 1);
            used_clean = false(C, 1);
            cMat = corrMat;

            for k = 1:C
                cMat(used_raw, :)   = -1;
                cMat(:, used_clean) = -1;
                [~, idx] = max(cMat(:));
                [ri, ci] = ind2sub([C C], idx);
                pairs(k,:)     = [ri, ci];
                used_raw(ri)   = true;
                used_clean(ci) = true;
            end
            [~, ord] = sort(pairs(:,1));
            pairs = pairs(ord, :);
        end

    end % static private methods

end
