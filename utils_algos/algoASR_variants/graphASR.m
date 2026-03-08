classdef graphASR < handle
    % graphASR - Graph Signal Processing Artifact Subspace Reconstruction.
    %   Detects and suppresses EEG artifacts by operating in the Graph Fourier
    %   domain. Channel topology is encoded as a k-NN correlation graph;
    %   artifact components are identified by inflated graph-frequency variance
    %   and suppressed via configurable shrinkage operators.
    %
    %   Constructor signature (matches Experimenting.m harness):
    %       obj = graphASR(srate, params)
    %       obj = graphASR(srate, params, tProbe)

    properties (Access = public)

        %% Core
        nchans
        srate
        nsamples = 0
        modifiedMask

        %% Graph Model
        W                   % Weighted adjacency matrix
        L                   % Graph Laplacian
        U                   % Graph Fourier basis (eigenvectors of L)
        lam                 % Graph frequencies (eigenvalues of L)
        Vp                  % Graph-PCA principal components
        r                   % Number of retained graph-PCA components
        vg0                 % Baseline per-frequency variance (from calibration)

        %% Parameters
        shrinkStyle = 'zero'    % Shrinkage mode: 'zero' | 'cap' | 'gamma'
        inflation   = 5         % Variance inflation threshold
        zeta        = 10        % Shrinkage aggressiveness (cap/gamma modes)
        varKeep     = 0.9999    % Fraction of graph variance retained in PCA

        %% Graph tuning
        kNeighbors  = 4         % k-NN degree for adjacency construction
        shrinkCorr  = 0.05      % Correlation shrinkage toward identity
        lambdaGain  = 5         % Frequency-dependent detection gain
        lambdaExp   = 0.5       % Exponent for frequency weighting curve

        %% Windowing
        window_len              % Window length in samples (set in constructor)
        hop_len                 % Hop length in samples (set in constructor)

        %% State
        last_trivial  = true    % Trivial flag from last processed window
        % FIX #4: persist first_iteration across process() calls to avoid
        % hard-write discontinuity at segment boundaries (cal → closed → open)
        is_first_window = true

        %% Telemetry
        tProbe
        probeRaw
        probeClean

    end


    methods

    %% ================= CONSTRUCTOR =================
    function obj = graphASR(srate, params, tProbe)
        % Constructor
        %   srate  : Sampling rate in Hz
        %   params : Struct of algorithm parameters (fields mapped to properties)
        %   tProbe : (optional) timeProbe instance

        if nargin < 1 || isempty(srate) || ~isnumeric(srate) || srate <= 0
            error('graphASR: srate must be a positive numeric scalar.');
        end

        obj.srate = srate;

        % FIX #9: window_len and hop_len set from srate here so that if
        % params overrides them, calibrate() will use the overridden values.
        obj.window_len = round(0.5  * srate);
        obj.hop_len    = round(0.25 * srate);

        if nargin >= 3 && ~isempty(tProbe)
            obj.tProbe = tProbe;
        else
            obj.tProbe = timeProbe();
        end

        if nargin >= 2 && isstruct(params)
            f = fieldnames(params);
            for i = 1:length(f)
                if isprop(obj, f{i})
                    obj.(f{i}) = params.(f{i});
                end
            end
        end

        obj.modifiedMask  = false(1, 0);
        obj.is_first_window = true;
    end


    %% ================= CALIBRATION =================
    function calibrate(obj, X)

        % FIX #8: input validation
        if isempty(X) || ~isnumeric(X)
            error('graphASR:calibrate', 'Calibration data must be a non-empty numeric matrix.');
        end

        fprintf('--- [graphASR] Graph calibration ---\n')
        obj.tProbe.start('calibration');

        [obj.nchans, ntime] = size(X);

        clip = 0.99999;

        % FIX #9: use obj.window_len / obj.hop_len instead of hardcoded values
        % so params overrides are respected during graph construction.
        wLen = obj.window_len;
        hLen = round(0.05 * obj.srate);     % fine-grained hop for correlation averaging

        starts = 1:hLen:(ntime - wLen + 1);

        if isempty(starts)
            error('graphASR:calibrate', ...
                'Calibration data too short for window_len=%d samples.', wLen);
        end

        zSum   = zeros(obj.nchans);
        zCount = zeros(obj.nchans);

        for t = starts
            Xw   = X(:, t:t+wLen-1);
            rho  = corr(Xw');
            rho  = max(min(rho, clip), -clip);
            mask = isfinite(rho);
            zSum(mask)   = zSum(mask)   + atanh(rho(mask));
            zCount(mask) = zCount(mask) + 1;
        end

        corrMat = tanh(zSum ./ max(zCount, 1));

        % Shrink toward identity, enforce symmetry, fix diagonal
        corrMat = (1 - obj.shrinkCorr) * corrMat + obj.shrinkCorr * eye(obj.nchans);
        corrMat = 0.5 * (corrMat + corrMat');
        corrMat(logical(eye(obj.nchans))) = 1;

        %% Graph adjacency (k-NN on absolute correlation)
        S = abs(corrMat);
        S(logical(eye(obj.nchans))) = 0;

        [~, order] = sort(S, 2, 'descend');
        adjMask    = false(obj.nchans);
        k          = min(obj.kNeighbors, obj.nchans - 1);

        for i = 1:obj.nchans
            adjMask(i, order(i, 1:k)) = true;
        end

        obj.W = 0.5 * ((S .* adjMask) + (S .* adjMask)');
        obj.L = diag(sum(obj.W, 2)) - obj.W;

        %% Graph Fourier transform
        [obj.U, obj.lam] = eig(obj.L, 'vector');
        [obj.lam, idx]   = sort(obj.lam, 'ascend');
        obj.U            = obj.U(:, idx);

        %% Graph PCA
        Xg       = obj.U' * X;
        graphCov = (Xg * Xg') / ntime;

        % FIX #10: enforce symmetry before eig to guard against floating-point asymmetry
        graphCov     = 0.5 * (graphCov + graphCov');
        [Vp, Dp]     = eig(graphCov, 'vector');
        [Dp, ord]    = sort(Dp, 'descend');
        obj.Vp       = Vp(:, ord);

        obj.r = find(cumsum(Dp) / sum(Dp) >= obj.varKeep, 1);
        if isempty(obj.r)
            obj.r = obj.nchans;
        end

        obj.vg0 = max(var(Xg, 0, 2), 1e-8);

        %% Probes
        obj.probeRaw   = graphASRprobe(12);
        obj.probeClean = graphASRprobe(12);
        obj.probeRaw.setBaseline(obj);
        obj.probeClean.setBaseline(obj);

        % Reset runtime state on re-calibration
        obj.nsamples      = 0;
        obj.modifiedMask  = false(1, 0);
        obj.last_trivial  = true;
        obj.is_first_window = true;

        obj.tProbe.stop('calibration');
    end


    %% ================= PROCESS =================
    function Xout = process(obj, X)

        % FIX #8: input validation
        if isempty(obj.U)
            error('graphASR:notCalibrated', 'Run calibrate() before process().');
        end
        if isempty(X) || ~isnumeric(X)
            error('graphASR:invalidInput', 'Data must be a non-empty numeric matrix.');
        end
        if size(X, 1) ~= obj.nchans
            error('graphASR:channelMismatch', ...
                'Expected %d channels, got %d.', obj.nchans, size(X,1));
        end

        obj.tProbe.start('process');

        [~, T] = size(X);

        % Pre-extend modifiedMask for this block
        obj.modifiedMask(obj.nsamples + (1:T)) = false;

        Xout     = X;
        Vp_clean = obj.Vp(:, 1:obj.r);
        P_graph  = Vp_clean * Vp_clean';

        lamNorm = (obj.lam - min(obj.lam)) ./ (max(obj.lam) - min(obj.lam) + eps);
        alpha   = 1 + obj.lambdaGain * (lamNorm .^ obj.lambdaExp);

        overlap = obj.window_len - obj.hop_len;

        for t = 1:obj.hop_len:(T - obj.window_len + 1)

            idx    = t : t + obj.window_len - 1;
            Xw_raw = X(:, idx);

            %% Graph transform & detection
            gSig          = obj.U' * Xw_raw;
            lambda_sorted = max(var(gSig, 0, 2), 1e-12);
            riem_dist     = norm(log(lambda_sorted ./ obj.vg0));

            rho_eff        = (lambda_sorted ./ obj.vg0) .* alpha;
            hi             = rho_eff > obj.inflation;
            trivial_window = ~any(hi);

            %% Shrinkage
            shrink = ones(size(rho_eff));
            if ~trivial_window
                switch lower(obj.shrinkStyle)
                    case 'zero'
                        shrink(hi) = 0;
                    case 'cap'
                        shrink(hi) = sqrt(obj.inflation ./ (obj.zeta * rho_eff(hi)));
                    case 'gamma'
                        shrink(hi) = (obj.inflation ./ rho_eff(hi)) .^ obj.zeta;
                end
            end

            %% Reconstruction in graph domain
            gSig_s = shrink .* gSig;
            if ~trivial_window
                gSig_s = P_graph * gSig_s;
            end
            Xw_cleaned = real(obj.U * gSig_s);

            %% Overlap-add blend
            % FIX #1/#4: is_first_window persists across process() calls.
            % On the first window ever (or first after calibrate()), do a
            % direct write. On all subsequent windows, blend the overlap zone
            % and direct-write the non-overlap zone.
            if obj.is_first_window
                Xout(:, idx)    = Xw_cleaned;
                obj.is_first_window = false;
            else
                % --- Overlap zone: raised-cosine blend ---
                blend     = (1 - cos(pi * (1:overlap) / overlap)) / 2;
                sub_range = idx(1:overlap);

                % FIX #3: only mark modified when reconstruction actually occurred
                if ~trivial_window || ~obj.last_trivial
                    obj.modifiedMask(obj.nsamples + sub_range) = true;
                end

                Xout(:, sub_range) = ...
                    blend       .* Xw_cleaned(:, 1:overlap) + ...
                    (1 - blend) .* Xout(:, sub_range);

                % --- Non-overlap zone: direct write ---
                % FIX #1: this was already a direct write in the original;
                % kept as-is because the non-overlap samples haven't been
                % touched by any prior window — they carry no prior estimate
                % to blend against. Mark as modified if reconstruction occurred.
                non_ol = idx(overlap+1:end);
                if ~trivial_window
                    obj.modifiedMask(obj.nsamples + non_ol) = true;
                end
                Xout(:, non_ol) = Xw_cleaned(:, overlap+1:end);
            end

            obj.last_trivial = trivial_window;

            %% RAW probe (always, every window)
            obj.probeRaw.update( ...
                lambda_sorted,          ...
                obj.vg0 * obj.inflation, ...
                [],                     ...
                trivial_window,         ...
                riem_dist);

            %% CLEAN probe
            % FIX #2: read from Xout AFTER the blend is complete so the
            % clean probe always measures a consistent, fully-written signal.
            % Only log when reconstruction occurred (matching emaASR convention).
            if ~trivial_window || ~obj.last_trivial
                gSig_clean   = obj.U' * Xout(:, idx);
                lambda_clean = max(var(gSig_clean, 0, 2), 1e-12);
                riem_clean   = norm(log(lambda_clean ./ obj.vg0));
                obj.probeClean.update( ...
                    lambda_clean,            ...
                    obj.vg0 * obj.inflation,  ...
                    [],                      ...
                    trivial_window,          ...
                    riem_clean);
            end

        end

        obj.nsamples = obj.nsamples + T;
        obj.tProbe.stop('process');

    end

    end % end methods
end
