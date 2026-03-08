classdef hmoASRprobe < asrProbe
    % hmoASRprobe - Diagnostic probe for hmoASR.
    %   Extends asrProbe with metrics specific to the HMO adaptive model:
    %     - Geodesic drift of the live clean covariance from calibration
    %     - Per-window rejection count and per-PC inflation
    %     - Threshold adaptation velocity relative to calibration anchor
    %     - Basis coherence (Modal Assurance Criterion diagonal)
    %     - Effective clean-sample accumulation
    %     - Step-to-step covariance drift
    %
    %   Shared metrics (from asrProbe):
    %     traceCov, maxEigCov, condCov, spectralEntropy, logDetDrift,
    %     riemannDrift, eigInflationMax/TopK, dirInflation*, energyRatio,
    %     normR, rankR, trivialFlag
    %
    %   Usage:
    %       p = hmoASRprobe(nchans);
    %       p.setBaseline(hmoObj);      % once, after calibrate()
    %       % inside each process window:
    %       p.update(hmoObj, Ct, Ek, lam, Tproj, R, reject, riem_dist);

    properties (Access = public)

        % ===== HMO Calibration Baseline (frozen at calibration) =====
        C_init          % C_ch at t=0
        V_init          % E   at t=0
        Tpc_init        % Tpc at t=0  [C x 1]

        % ===== Model Adaptation Tracking =====
        modelRiemannDrift   % AIRM dist(C_init, C_ch_t) — clean model drift
        thresholdAdaptation % ||Tpc_t - Tpc_0|| / ||Tpc_0|| — relative threshold shift
        sampleWeightTotal   % n_tot at each window

        % ===== Basis Coherence =====
        basisCoherence      % Mean diagonal of MAC(V_init, E_t)

        % ===== Block-wise Detection =====
        riemannDistRaw      % Geodesic dist(C_ch_t, C_block) per window
        rejectionCount      % Number of PCs rejected per window
        maxInflation        % max(lambda / Tproj) per window
        meanInflation       % mean(lambda / Tproj) per window

        % ===== Step-to-step Drift =====
        riemannDriftStep    % ||log(lambda_t ./ lambda_{t-1})||

    end

    properties (Access = private)
        prevLambda_hmo      % For step-to-step drift
    end

    % ================================================================
    methods

        function obj = hmoASRprobe(subspaceRank)
            if nargin < 1, subspaceRank = []; end
            obj = obj@asrProbe(subspaceRank);
            obj.resetOwn();
        end

        % ── Baseline ─────────────────────────────────────────────────
        function setBaseline(obj, hmoObj)
            % SETBASELINE - Freeze calibration state as the drift reference.
            %   hmoObj : hmoASR (or hebbASR) instance after calibrate()
            if ~obj.enabled, return; end

            if ~isprop(hmoObj, 'C_ch') || isempty(hmoObj.C_ch)
                error('hmoASRprobe:setBaseline', ...
                    'hmoObj must have a non-empty property C_ch.');
            end

            obj.C_init   = hmoObj.C_ch;
            obj.V_init   = hmoObj.E;
            obj.Tpc_init = hmoObj.Tpc;

            if isempty(obj.subspaceRank)
                obj.subspaceRank = size(hmoObj.C_ch, 1);
            end

            % Populate base-class lambda0, V0, C0 from C_ch
            C0_sym = (obj.C_init + obj.C_init.') / 2;
            [V0, D0] = eig(C0_sym);
            lambda   = diag(D0);  lambda(lambda < 1e-12) = 1e-12;
            [lam_sorted, idx] = sort(lambda, 'descend');
            V0 = V0(:, idx);

            obj.lambda0     = lam_sorted(:);
            obj.V0          = V0(:, 1:min(obj.subspaceRank, numel(lam_sorted)));
            obj.C0          = C0_sym;
            obj.C0_inv_sqrt = V0 * diag(1./sqrt(lam_sorted)) * V0.';

            if ~isempty(hmoObj.Tpc)
                obj.T0 = hmoObj.Tpc(:);
            end
        end

        % ── Update ───────────────────────────────────────────────────
        function update(obj, hmoObj, Ct, Ek, lam, Tproj, R, reject, riem_dist)
            % UPDATE - Record per-window metrics.
            %   hmoObj    : current hmoASR/hebbASR instance
            %   Ct        : [C x C] current block covariance
            %   Ek        : [C x C] current eigenbasis
            %   lam       : [C x 1] current eigenvalues
            %   Tproj     : [C x 1] projected per-PC thresholds
            %   R         : [C x C] reconstruction (projector) matrix
            %   reject    : [C x 1] logical rejection mask
            %   riem_dist : scalar AIRM distance (block vs clean model)

            if ~obj.enabled, return; end

            lam_col  = double(lam(:));
            lam_col(lam_col < 1e-12) = 1e-12;
            lam_desc = sort(lam_col, 'descend');

            % 1. Base-class shared metrics
            update@asrProbe(obj, lam_desc, Tproj, R, ~any(reject), riem_dist);

            % 2. Model geodesic drift from calibration anchor
            obj.modelRiemannDrift(end+1,1) = ...
                obj.calcGeodesicSPD(obj.C_init, hmoObj.C_ch);

            % 3. Effective sample count
            obj.sampleWeightTotal(end+1,1) = hmoObj.n_tot;

            % 4. Threshold adaptation velocity
            if ~isempty(obj.Tpc_init) && ~isempty(hmoObj.Tpc)
                obj.thresholdAdaptation(end+1,1) = ...
                    norm(hmoObj.Tpc(:) - obj.Tpc_init(:)) / ...
                    (norm(obj.Tpc_init(:)) + eps);
            else
                obj.thresholdAdaptation(end+1,1) = NaN;
            end

            % 5. Basis coherence — Modal Assurance Criterion diagonal mean
            if ~isempty(obj.V_init) && ~isempty(Ek)
                k   = min(size(obj.V_init,2), size(Ek,2));
                MAC = (obj.V_init(:,1:k).' * Ek(:,1:k)) .^ 2;
                obj.basisCoherence(end+1,1) = mean(diag(MAC));
            else
                obj.basisCoherence(end+1,1) = NaN;
            end

            % 6. Block-wise detection stats
            obj.riemannDistRaw(end+1,1) = riem_dist;
            obj.rejectionCount(end+1,1) = sum(reject);

            Tproj_col = double(Tproj(:));
            n_min     = min(numel(lam_col), numel(Tproj_col));
            inflation = lam_col(1:n_min) ./ max(Tproj_col(1:n_min), eps);
            obj.maxInflation(end+1,1)  = max(inflation);
            obj.meanInflation(end+1,1) = mean(inflation);

            % 7. Step-to-step covariance drift
            if ~isempty(obj.prevLambda_hmo)
                k = min(numel(obj.prevLambda_hmo), numel(lam_desc));
                obj.riemannDriftStep(end+1,1) = ...
                    norm(log(lam_desc(1:k) ./ obj.prevLambda_hmo(1:k)));
            else
                obj.riemannDriftStep(end+1,1) = NaN;
            end
            obj.prevLambda_hmo = lam_desc;
        end

        % ── Reset ────────────────────────────────────────────────────
        function reset(obj)
            reset@asrProbe(obj);
            obj.resetOwn();
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function resetOwn(obj)
            obj.prevLambda_hmo      = [];
            obj.modelRiemannDrift   = zeros(0,1);
            obj.thresholdAdaptation = zeros(0,1);
            obj.sampleWeightTotal   = zeros(0,1);
            obj.basisCoherence      = zeros(0,1);
            obj.riemannDistRaw      = zeros(0,1);
            obj.rejectionCount      = zeros(0,1);
            obj.maxInflation        = zeros(0,1);
            obj.meanInflation       = zeros(0,1);
            obj.riemannDriftStep    = zeros(0,1);
        end

        function dist = calcGeodesicSPD(~, C1, C2)
            % CALCGEODESICSPD - Affine-Invariant Riemannian Metric.
            if isempty(C1) || isempty(C2), dist = NaN; return; end
            try
                C1s    = (C1 + C1.') / 2;
                [V, D] = eig(C1s);
                d      = diag(D);  d(d < eps) = eps;
                i_sqrt = V * diag(1./sqrt(d)) * V.';
                mat    = i_sqrt * C2 * i_sqrt;
                vals   = eig((mat + mat.')/2);
                dist   = sqrt(sum(log(max(eps, vals)).^2));
            catch
                dist = NaN;
            end
        end

    end % private methods

end
