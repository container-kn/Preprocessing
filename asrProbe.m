classdef asrProbe < handle
    % asrProbe - Base class for ASR diagnostic telemetry.
    %   Captures per-window metrics that are meaningful across ALL ASR variants:
    %   covariance structure, spectral entropy, inflation, reconstruction quality,
    %   and Riemannian drift from a calibration baseline.
    %
    %   Subclasses extend update() with variant-specific arguments:
    %       vanillaASRprobe  < asrProbe  — adds subspace angles, threshold drift
    %       graphASRprobe    < asrProbe  — adds graph-spectral metrics, step drift
    %       hmoASRprobe      < asrProbe  — adds model adaptation tracking
    %
    %   Usage (base):
    %       p = asrProbe(12);
    %       p.setBaseline(asrObj);          % captures lambda0 from M*M'
    %       p.update(lambda, thr, R, trivial, riem_dist);
    %       p.reset();

    properties (Access = public)
        enabled      logical = true     % Global toggle — false = zero overhead
        step    (1,1) double = 0        % Window counter
        subspaceRank                    % Number of principal components tracked

        % ===== Calibration Baseline =====
        lambda0         % Baseline eigenvalues (descending) from M*M'
        C0              % Baseline covariance  M*M'
        C0_inv_sqrt     % C0^{-1/2} for AIRM computation
        V0              % Baseline eigenvectors (top subspaceRank columns)
        T0              % Baseline threshold diagonal (from state.T)

        % ===== Covariance Metrics =====
        traceCov        % trace(Ct)
        maxEigCov       % lambda_1
        condCov         % lambda_1 / lambda_C

        % ===== Spectral Structure =====
        spectralEntropy % Shannon entropy of normalised eigenspectrum
        logDetDrift     % sum(log lambda_t) - sum(log lambda_0)

        % ===== Riemannian Distance =====
        riemannDrift    % AIRM distance from calibration baseline

        % ===== Global Eigen Inflation =====
        eigInflationMax     % max(lambda_t / lambda_0)
        eigInflationTopK    % mean inflation over top-k components

        % ===== Directional Inflation =====
        dirInflationMax     % max(lambda / threshold_dir)
        dirInflationMean    % mean(lambda / threshold_dir)
        dirNumViolations    % sum(lambda > threshold_dir)

        % ===== Energy Ratio =====
        energyRatio     % Fraction of total energy in baseline subspace

        % ===== Reconstruction =====
        normR           % ||R - I||_F
        rankR           % rank(R)
        trivialFlag     % 1 = clean window, 0 = reconstruction occurred
    end

    methods

        % =====================================================
        % Constructor
        % =====================================================
        function obj = asrProbe(subspaceRank)
            if nargin < 1, subspaceRank = []; end
            obj.subspaceRank = subspaceRank;
            obj.reset();
        end

        % =====================================================
        % Baseline Capture
        % =====================================================
        function setBaseline(obj, state0)
            % SETBASELINE - Capture calibration state as the drift anchor.
            %   state0 must be an object with property M (mixing matrix).
            %   Optionally reads state0.T for threshold drift tracking.
            if ~obj.enabled, return; end

            if ~isprop(state0, 'M') || isempty(state0.M)
                error('asrProbe:setBaseline', ...
                    'state0 must have a non-empty property M.');
            end

            obj.C0 = (state0.M * state0.M' + state0.M * state0.M') / 2; % enforce symmetry

            [V, D]   = eig(obj.C0);
            lambda   = diag(D);
            lambda(lambda < 1e-12) = 1e-12;

            [lambda_sorted, idx] = sort(lambda, 'descend');
            V = V(:, idx);

            obj.lambda0 = lambda_sorted(:);

            if isempty(obj.subspaceRank)
                obj.subspaceRank = length(lambda_sorted);
            end

            k        = min(obj.subspaceRank, length(lambda_sorted));
            obj.V0   = V(:, 1:k);

            % C0^{-1/2} for AIRM — use full V, not truncated
            obj.C0_inv_sqrt = V * diag(1 ./ sqrt(lambda_sorted)) * V';

            % Threshold baseline — diag of T operator
            if isprop(state0, 'T') && ~isempty(state0.T)
                obj.T0 = diag(state0.T);
            elseif isfield(state0, 'T') && ~isempty(state0.T)
                obj.T0 = diag(state0.T);
            else
                obj.T0 = [];
            end
        end

        % =====================================================
        % Core Update — shared metrics only
        % =====================================================
        function update(obj, lambda, threshold_dir, R, trivial, riem_dist)
            % UPDATE - Record per-window metrics shared across all variants.
            %   lambda        : [C x 1] eigenvalues of current window covariance
            %   threshold_dir : [C x 1] per-component thresholds (or [])
            %   R             : [C x C] reconstruction matrix (or [])
            %   trivial       : logical scalar — true if R == I
            %   riem_dist     : scalar AIRM distance from baseline
            %
            %   Subclasses call: update@asrProbe(obj,...) then add their own.

            if ~obj.enabled, return; end

            obj.step = obj.step + 1;

            lambda = obj.sanitiseLambda(lambda);

            % --- Covariance metrics ---
            obj.traceCov(end+1,1)  = sum(lambda);
            obj.maxEigCov(end+1,1) = lambda(1);
            obj.condCov(end+1,1)   = lambda(1) / lambda(end);

            % --- Spectral entropy ---
            p = lambda / sum(lambda);
            p(p < 1e-12) = 1e-12;
            obj.spectralEntropy(end+1,1) = -sum(p .* log(p));

            % --- Log-det drift ---
            if ~isempty(obj.lambda0)
                k = min(length(lambda), length(obj.lambda0));
                obj.logDetDrift(end+1,1) = ...
                    sum(log(lambda(1:k))) - sum(log(obj.lambda0(1:k)));
            else
                obj.logDetDrift(end+1,1) = NaN;
            end

            % --- Riemannian drift ---
            obj.riemannDrift(end+1,1) = riem_dist;

            % --- Global eigen inflation ---
            if ~isempty(obj.lambda0)
                k          = min(length(lambda), length(obj.lambda0));
                inflation  = lambda(1:k) ./ obj.lambda0(1:k);
                obj.eigInflationMax(end+1,1)  = max(inflation);
                kk = min(obj.subspaceRank, k);
                obj.eigInflationTopK(end+1,1) = mean(inflation(1:kk));
            else
                obj.eigInflationMax(end+1,1)  = NaN;
                obj.eigInflationTopK(end+1,1) = NaN;
            end

            % --- Directional inflation ---
            if ~isempty(threshold_dir)
                threshold_dir = threshold_dir(:);
                n_min = min(length(lambda), length(threshold_dir));
                inflation_dir = lambda(1:n_min) ./ threshold_dir(1:n_min);
                obj.dirInflationMax(end+1,1)  = max(inflation_dir);
                obj.dirInflationMean(end+1,1) = mean(inflation_dir);
                obj.dirNumViolations(end+1,1) = sum(inflation_dir > 1);
            else
                obj.dirInflationMax(end+1,1)  = NaN;
                obj.dirInflationMean(end+1,1) = NaN;
                obj.dirNumViolations(end+1,1) = NaN;
            end

            % --- Energy ratio ---
            if ~isempty(obj.lambda0)
                k = min(obj.subspaceRank, length(lambda));
                obj.energyRatio(end+1,1) = sum(lambda(1:k)) / sum(lambda);
            else
                obj.energyRatio(end+1,1) = NaN;
            end

            % --- Reconstruction ---
            if ~isempty(R)
                obj.normR(end+1,1) = norm(R - eye(size(R)), 'fro');
                obj.rankR(end+1,1) = rank(R);
            else
                obj.normR(end+1,1) = NaN;
                obj.rankR(end+1,1) = NaN;
            end

            obj.trivialFlag(end+1,1) = trivial;
        end

        % =====================================================
        % Reset
        % =====================================================
        function reset(obj)
            % RESET - Clear all accumulated metrics. Preserves baseline.
            obj.step = 0;

            obj.traceCov         = zeros(0,1);
            obj.maxEigCov        = zeros(0,1);
            obj.condCov          = zeros(0,1);
            obj.spectralEntropy  = zeros(0,1);
            obj.logDetDrift      = zeros(0,1);
            obj.riemannDrift     = zeros(0,1);
            obj.eigInflationMax  = zeros(0,1);
            obj.eigInflationTopK = zeros(0,1);
            obj.dirInflationMax  = zeros(0,1);
            obj.dirInflationMean = zeros(0,1);
            obj.dirNumViolations = zeros(0,1);
            obj.energyRatio      = zeros(0,1);
            obj.normR            = zeros(0,1);
            obj.rankR            = zeros(0,1);
            obj.trivialFlag      = zeros(0,1);
        end

    end

    methods (Access = protected)

        function lambda = sanitiseLambda(~, lambda)
            % SANITISELAMBDA - Enforce column vector, floor at 1e-12, sort descending.
            lambda = double(lambda(:));
            lambda(lambda < 1e-12) = 1e-12;
            lambda = sort(lambda, 'descend');
        end

    end

end
