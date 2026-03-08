classdef hebbASRprobe < hmoASRprobe
    % hebbASRprobe - Diagnostic probe for hebbASR (PSP/PSW).
    %   Extends hmoASRprobe with metrics specific to Hebbian weight learning:
    %
    %   Hebbian-specific metrics:
    %     weightDeltaW    — ||W_t - W_{t-1}||_F / ||W_{t-1}||_F
    %                       Feedforward weight velocity
    %     weightDeltaM    — ||M_t - M_{t-1}||_F / ||M_{t-1}||_F
    %                       Lateral weight velocity
    %     weightCovAlign  — ||W'W - C_ch||_F / ||C_ch||_F
    %                       Convergence quality: W should track sqrt(C_ch) at convergence
    %     mCondition      — cond(M) — tracks near-singularity of lateral weights
    %     hebbStepCount   — obj.hebb_step value at each window
    %
    %   Inherited from hmoASRprobe:
    %     modelRiemannDrift, thresholdAdaptation, sampleWeightTotal,
    %     basisCoherence, riemannDistRaw, rejectionCount, maxInflation,
    %     meanInflation, riemannDriftStep
    %
    %   Inherited from asrProbe (base):
    %     traceCov, maxEigCov, condCov, spectralEntropy, logDetDrift,
    %     riemannDrift, eigInflation*, dirInflation*, energyRatio,
    %     normR, rankR, trivialFlag
    %
    %   Usage:
    %       p = hebbASRprobe(nchans);
    %       p.setBaseline(hebbObj);     % once, after calibrate()
    %       p.update(hebbObj, Ct, Ek, lam, Tproj, R, reject, riem_dist);

    properties (Access = public)

        % ===== Hebbian Baseline (frozen at calibration) =====
        W_init      % W at t=0  [C x C]
        M_init      % M at t=0  [C x C]

        % ===== Feedforward Weight Velocity =====
        weightDeltaW    % ||W_t - W_{t-1}||_F / ||W_{t-1}||_F

        % ===== Lateral Weight Velocity =====
        weightDeltaM    % ||M_t - M_{t-1}||_F / ||M_{t-1}||_F

        % ===== Weight-Covariance Alignment =====
        weightCovAlign  % ||W'W - C_ch||_F / ||C_ch||_F

        % ===== Lateral Matrix Conditioning =====
        mCondition      % cond(M)

        % ===== Learning Dynamics =====
        hebbStepCount   % hebb_step at each window

    end

    properties (Access = private)
        prevW       % W from previous window
        prevM       % M from previous window
    end

    % ================================================================
    methods

        function obj = hebbASRprobe(subspaceRank)
            if nargin < 1, subspaceRank = []; end
            obj = obj@hmoASRprobe(subspaceRank);
            obj.resetOwn();
        end

        % ── Baseline ─────────────────────────────────────────────────
        function setBaseline(obj, hebbObj)
            % SETBASELINE - Freeze calibration state including Hebbian weights.
            %   Calls hmoASRprobe.setBaseline for all HMO metrics, then captures W, M.

            setBaseline@hmoASRprobe(obj, hebbObj);

            if isprop(hebbObj, 'W') && ~isempty(hebbObj.W)
                obj.W_init  = hebbObj.W;
                obj.prevW   = hebbObj.W;
            end
            if isprop(hebbObj, 'M') && ~isempty(hebbObj.M)
                obj.M_init  = hebbObj.M;
                obj.prevM   = hebbObj.M;
            end
        end

        % ── Update ───────────────────────────────────────────────────
        function update(obj, hebbObj, Ct, Ek, lam, Tproj, R, reject, riem_dist)
            % UPDATE - Record all metrics for hebbASR.
            %   hebbObj   : current hebbASR instance
            %   Ct        : [C x C] current block covariance
            %   Ek        : [C x C] current eigenbasis
            %   lam       : [C x 1] current eigenvalues
            %   Tproj     : [C x 1] projected per-PC thresholds
            %   R         : [C x C] reconstruction matrix
            %   reject    : [C x 1] logical rejection mask
            %   riem_dist : scalar AIRM distance from calibration

            if ~obj.enabled, return; end

            % 1. Parent HMO metrics (handles asrProbe base class chain too)
            update@hmoASRprobe(obj, hebbObj, Ct, Ek, lam, Tproj, R, reject, riem_dist);

            % 2. Current weight matrices
            W_cur = hebbObj.W;
            M_cur = hebbObj.M;

            % 3. Feedforward weight velocity
            if ~isempty(obj.prevW) && ~isempty(W_cur)
                obj.weightDeltaW(end+1,1) = ...
                    norm(W_cur - obj.prevW, 'fro') / max(norm(obj.prevW, 'fro'), eps);
            else
                obj.weightDeltaW(end+1,1) = NaN;
            end

            % 4. Lateral weight velocity
            if ~isempty(obj.prevM) && ~isempty(M_cur)
                obj.weightDeltaM(end+1,1) = ...
                    norm(M_cur - obj.prevM, 'fro') / max(norm(obj.prevM, 'fro'), eps);
            else
                obj.weightDeltaM(end+1,1) = NaN;
            end

            obj.prevW = W_cur;
            obj.prevM = M_cur;

            % 5. Weight-covariance alignment
            %    At convergence W'W → C_ch (PSP) or W → sqrt(C_ch)·E (PSW).
            %    Track W'W vs C_ch as a convergence diagnostic.
            if ~isempty(W_cur) && ~isempty(hebbObj.C_ch)
                WtW    = W_cur.' * W_cur;
                C_norm = norm(hebbObj.C_ch, 'fro');
                obj.weightCovAlign(end+1,1) = ...
                    norm(WtW - hebbObj.C_ch, 'fro') / max(C_norm, eps);
            else
                obj.weightCovAlign(end+1,1) = NaN;
            end

            % 6. Lateral matrix conditioning
            if ~isempty(M_cur)
                try
                    obj.mCondition(end+1,1) = cond(M_cur);
                catch
                    obj.mCondition(end+1,1) = NaN;
                end
            else
                obj.mCondition(end+1,1) = NaN;
            end

            % 7. Hebbian step counter
            if isprop(hebbObj, 'hebb_step')
                obj.hebbStepCount(end+1,1) = hebbObj.hebb_step;
            else
                obj.hebbStepCount(end+1,1) = NaN;
            end
        end

        % ── Reset ────────────────────────────────────────────────────
        function reset(obj)
            reset@hmoASRprobe(obj);
            obj.resetOwn();
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function resetOwn(obj)
            obj.prevW          = [];
            obj.prevM          = [];
            obj.weightDeltaW   = zeros(0,1);
            obj.weightDeltaM   = zeros(0,1);
            obj.weightCovAlign = zeros(0,1);
            obj.mCondition     = zeros(0,1);
            obj.hebbStepCount  = zeros(0,1);
        end

    end % private methods

end
