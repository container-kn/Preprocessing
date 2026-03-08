classdef emaASRprobe < handle
    % emaASRprobe - Diagnostic telemetry for emaASR.
    %   Extends the asrProbe metric set with EMA-specific adaptation metrics:
    %   manifold velocity (adaptationDeltaM/T), buffer fill fraction, and
    %   stepwise subspace angle — all measured against the frozen calibration
    %   anchor so drift is interpretable in absolute terms.
    %
    %   Metric groups (aligned with asrProbe convention):
    %     Baseline          : C0, V0, lambda0, T0 (set once via setBaseline)
    %     Covariance        : traceCov, maxEigCov, condCov
    %     Spectral          : spectralEntropy, logDetDrift
    %     Drift             : riemannDrift, riemannDriftStep, energyRatio,
    %                         thresholdDrift, logDetDrift
    %     Inflation         : eigInflationMax, eigInflationTopK,
    %                         dirInflationMax, dirInflationMean, dirNumViolations
    %     Subspace geometry : subspaceAngleMax, subspaceAngleMean,
    %                         subspaceAngleStep
    %     Reconstruction    : normR, rankR, trivialFlag
    %     EMA-specific      : adaptationDeltaM, adaptationDeltaT,
    %                         cleanBufferFill

    properties
        enabled      logical = true
        step         (1,1) double = 0
        subspaceRank

        % ===== Baseline (anchor — frozen at calibration) =====
        C0            % M*M' at calibration
        C0_inv_sqrt   % For AIRM distance computation
        V0            % Principal components at calibration
        lambda0       % Eigenvalues at calibration (descending)
        T0            % Threshold diagonal at calibration

        % ===== Covariance =====
        traceCov
        maxEigCov
        condCov

        % ===== Spectral structure =====
        spectralEntropy
        logDetDrift

        % ===== Drift =====
        riemannDrift      % AIRM distance from calibration anchor
        riemannDriftStep  % Step-wise eigenvalue drift (log-ratio norm)
        energyRatio       % Fraction of total energy in calibration subspace
        thresholdDrift    % Relative change in T from calibration anchor

        % ===== Inflation =====
        eigInflationMax
        eigInflationTopK
        dirInflationMax
        dirInflationMean
        dirNumViolations

        % ===== Subspace geometry =====
        subspaceAngleMax   % Max principal angle between Vt and V0
        subspaceAngleMean  % Mean principal angle between Vt and V0
        subspaceAngleStep  % Step-wise angle: Vt vs V_{t-1}

        % ===== Reconstruction =====
        normR
        rankR
        trivialFlag

        % ===== EMA-specific =====
        adaptationDeltaM  % ||M_t - M_{t-1}||_F / ||M_{t-1}||_F — manifold velocity
        adaptationDeltaT  % ||T_t - T_{t-1}||_F / ||T_{t-1}||_F — threshold velocity
        cleanBufferFill   % Clean buffer fill fraction (0-1)

        % ===== Internal =====
        prevLambda        % Lambda from previous step for riemannDriftStep
        prevV             % V from previous step for subspaceAngleStep
    end

    methods

        % =====================================================
        % Constructor
        % =====================================================
        function obj = emaASRprobe(subspaceRank)
            if nargin < 1
                subspaceRank = [];
            end
            obj.subspaceRank = subspaceRank;
            obj.reset();
        end

        % =====================================================
        % Baseline — called once after calibrate()
        % =====================================================
        function setBaseline(obj, asrObj)
            % SETBASELINE - Freeze calibration state as the drift reference.
            %   asrObj : emaASR instance (must have .M and .T populated)

            if ~obj.enabled, return; end

            if ~isprop(asrObj, 'M') || isempty(asrObj.M)
                error('emaASRprobe:setBaseline', ...
                    'asrObj must have a populated .M property.');
            end

            obj.C0     = asrObj.M * asrObj.M';
            obj.C0     = (obj.C0 + obj.C0') / 2;

            [V, D]     = eig(obj.C0);
            lambda     = diag(D);
            lambda(lambda < 1e-12) = 1e-12;

            [lambda_sorted, idx] = sort(lambda, 'descend');
            V = V(:, idx);

            obj.lambda0 = lambda_sorted(:);

            if isempty(obj.subspaceRank)
                obj.subspaceRank = numel(lambda_sorted);
            end

            k          = min(obj.subspaceRank, numel(lambda_sorted));
            obj.V0     = V(:, 1:k);

            obj.C0_inv_sqrt = V * diag(1 ./ sqrt(lambda_sorted)) * V';

            if isprop(asrObj, 'T') && ~isempty(asrObj.T)
                obj.T0 = diag(asrObj.T);
            else
                obj.T0 = [];
            end
        end

        % =====================================================
        % Update — called each processing window
        % =====================================================
        function update(obj, Ct, Vc, lambda_sorted, R, trivial, ...
                        riem_dist, M_current, T_current, ...
                        M_prev, T_prev, bufferFillFrac)
            % UPDATE - Record telemetry for one processing window.
            %
            %   Ct             : [C x C] current covariance
            %   Vc             : [C x C] current eigenvectors (columns)
            %   lambda_sorted  : [C x 1] eigenvalues (ascending, from eig sort)
            %   R              : [C x C] reconstruction matrix
            %   trivial        : logical scalar
            %   riem_dist      : AIRM distance from calibration anchor (scalar)
            %   M_current      : current EMA mixing matrix
            %   T_current      : current EMA threshold matrix
            %   M_prev         : mixing matrix from previous step
            %   T_prev         : threshold matrix from previous step
            %   bufferFillFrac : clean buffer fill fraction (0-1)

            if ~obj.enabled || isempty(Ct), return; end

            % Default optional args
            if nargin < 8,  M_current      = [];  end
            if nargin < 9,  T_current      = [];  end
            if nargin < 10, M_prev         = [];  end
            if nargin < 11, T_prev         = [];  end
            if nargin < 12, bufferFillFrac = NaN; end

            obj.step = obj.step + 1;

            lambda_sorted = double(lambda_sorted(:));
            lambda_sorted(lambda_sorted < 1e-12) = 1e-12;
            % emaASR returns eigenvalues ascending — reverse for consistency
            % with asrProbe convention (descending = largest first)
            lambda_d = flipud(lambda_sorted);

            % ===== Covariance =====
            obj.traceCov(end+1,1)  = sum(lambda_d);
            obj.maxEigCov(end+1,1) = lambda_d(1);
            obj.condCov(end+1,1)   = lambda_d(1) / lambda_d(end);

            % ===== Spectral entropy =====
            p = lambda_d / sum(lambda_d);
            p(p < 1e-12) = 1e-12;
            obj.spectralEntropy(end+1,1) = -sum(p .* log(p));

            % ===== Log-det drift =====
            if ~isempty(obj.lambda0)
                k = min(numel(lambda_d), numel(obj.lambda0));
                obj.logDetDrift(end+1,1) = ...
                    sum(log(lambda_d(1:k))) - sum(log(obj.lambda0(1:k)));
            else
                obj.logDetDrift(end+1,1) = NaN;
            end

            % ===== Subspace angles vs anchor =====
            if ~isempty(obj.V0) && ~isempty(Vc)
                k  = min(obj.subspaceRank, size(Vc, 2));
                sv = svd(obj.V0' * Vc(:, 1:k));
                sv = min(max(sv, -1), 1);
                angles = acos(sv);
                obj.subspaceAngleMax(end+1,1)  = max(angles);
                obj.subspaceAngleMean(end+1,1) = mean(angles);
            else
                obj.subspaceAngleMax(end+1,1)  = NaN;
                obj.subspaceAngleMean(end+1,1) = NaN;
            end

            % ===== Subspace angle vs previous step =====
            if ~isempty(obj.prevV) && ~isempty(Vc)
                k  = min(size(obj.prevV, 2), size(Vc, 2));
                sv = svd(obj.prevV(:,1:k)' * Vc(:,1:k));
                sv = min(max(sv, -1), 1);
                obj.subspaceAngleStep(end+1,1) = acos(sv(1));
            else
                obj.subspaceAngleStep(end+1,1) = NaN;
            end
            obj.prevV = Vc;

            % ===== Energy ratio in calibration subspace =====
            if ~isempty(obj.V0)
                obj.energyRatio(end+1,1) = trace(obj.V0' * Ct * obj.V0) / trace(Ct);
            else
                obj.energyRatio(end+1,1) = NaN;
            end

            % ===== Riemann drift from anchor =====
            obj.riemannDrift(end+1,1) = riem_dist;

            % ===== Step-wise Riemann drift (eigenvalue log-ratio) =====
            if ~isempty(obj.prevLambda)
                k = min(numel(obj.prevLambda), numel(lambda_d));
                obj.riemannDriftStep(end+1,1) = ...
                    norm(log(lambda_d(1:k) ./ obj.prevLambda(1:k)));
            else
                obj.riemannDriftStep(end+1,1) = NaN;
            end
            obj.prevLambda = lambda_d;

            % ===== Threshold drift vs anchor =====
            if ~isempty(obj.T0) && ~isempty(T_current)
                Tt = diag(T_current);
                obj.thresholdDrift(end+1,1) = norm(Tt - obj.T0) / norm(obj.T0);
            else
                obj.thresholdDrift(end+1,1) = NaN;
            end

            % ===== Eigenvalue inflation vs anchor =====
            if ~isempty(obj.lambda0)
                k          = min(numel(lambda_d), numel(obj.lambda0));
                inflation  = lambda_d(1:k) ./ obj.lambda0(1:k);
                obj.eigInflationMax(end+1,1) = max(inflation);
                kk = min(obj.subspaceRank, k);
                obj.eigInflationTopK(end+1,1) = mean(inflation(1:kk));
            else
                obj.eigInflationMax(end+1,1)  = NaN;
                obj.eigInflationTopK(end+1,1) = NaN;
            end

            % ===== Directional inflation vs threshold =====
            % threshold_dir not passed from emaASR — record NaN
            obj.dirInflationMax(end+1,1)  = NaN;
            obj.dirInflationMean(end+1,1) = NaN;
            obj.dirNumViolations(end+1,1) = NaN;

            % ===== Reconstruction =====
            if ~isempty(R)
                obj.normR(end+1,1) = norm(R - eye(size(R)), 'fro');
                obj.rankR(end+1,1) = rank(R);
            else
                obj.normR(end+1,1) = NaN;
                obj.rankR(end+1,1) = NaN;
            end
            obj.trivialFlag(end+1,1) = trivial;

            % ===== EMA-specific: adaptation velocity =====
            if ~isempty(M_prev) && ~isempty(M_current)
                obj.adaptationDeltaM(end+1,1) = ...
                    norm(M_current - M_prev, 'fro') / max(norm(M_prev, 'fro'), eps);
            else
                obj.adaptationDeltaM(end+1,1) = NaN;
            end

            if ~isempty(T_prev) && ~isempty(T_current)
                obj.adaptationDeltaT(end+1,1) = ...
                    norm(T_current - T_prev, 'fro') / max(norm(T_prev, 'fro'), eps);
            else
                obj.adaptationDeltaT(end+1,1) = NaN;
            end

            % ===== EMA-specific: buffer fill =====
            obj.cleanBufferFill(end+1,1) = bufferFillFrac;
        end

        % =====================================================
        % Reset
        % =====================================================
        function reset(obj)
            obj.step       = 0;
            obj.prevLambda = [];
            obj.prevV      = [];

            obj.traceCov        = zeros(0,1);
            obj.maxEigCov       = zeros(0,1);
            obj.condCov         = zeros(0,1);

            obj.spectralEntropy = zeros(0,1);
            obj.logDetDrift     = zeros(0,1);

            obj.riemannDrift    = zeros(0,1);
            obj.riemannDriftStep= zeros(0,1);
            obj.energyRatio     = zeros(0,1);
            obj.thresholdDrift  = zeros(0,1);

            obj.eigInflationMax  = zeros(0,1);
            obj.eigInflationTopK = zeros(0,1);

            obj.dirInflationMax  = zeros(0,1);
            obj.dirInflationMean = zeros(0,1);
            obj.dirNumViolations = zeros(0,1);

            obj.subspaceAngleMax  = zeros(0,1);
            obj.subspaceAngleMean = zeros(0,1);
            obj.subspaceAngleStep = zeros(0,1);

            obj.normR       = zeros(0,1);
            obj.rankR       = zeros(0,1);
            obj.trivialFlag = zeros(0,1);

            obj.adaptationDeltaM = zeros(0,1);
            obj.adaptationDeltaT = zeros(0,1);
            obj.cleanBufferFill  = zeros(0,1);
        end

    end
end
