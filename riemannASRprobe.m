classdef riemannASRprobe < asrProbe
    % riemannASRprobe - Diagnostic telemetry for riemannASR.
    %   Extends asrProbe with metrics specific to the Riemannian covariance
    %   averaging pipeline:
    %     - Karcher mean drift: AIRM distance between raw SCM and Karcher-averaged Ct
    %     - Step-to-step covariance drift: AIRM(Ct, Ct_prev)
    %     - Nonlinear eigenspace alignment: principal angles between Vt and V0
    %     - Threshold margin: per-component distance to threshold
    %
    %   update() signature:
    %       p.update(Ct, SCM, lambda, threshold_dir, R, trivial, riem_dist)
    %
    %   Ct            : Karcher-averaged covariance at this window [C x C]
    %   SCM           : Raw sample covariance at this window [C x C]
    %   lambda        : Sorted eigenvalues of Ct [C x 1]
    %   threshold_dir : Per-component threshold (diag of T*V product) [C x 1]
    %   R             : Reconstruction matrix [C x C] ([] if trivial)
    %   trivial       : logical — true if no reconstruction occurred
    %   riem_dist     : Precomputed AIRM distance from calibration baseline

    properties (Access = public)

        % ===== Karcher Mean vs Raw SCM =====
        karcherVsScmDist    % AIRM(Ct, SCM) — how much geodesic smoothing shifted the estimate

        % ===== Step-to-step drift =====
        riemannDriftStep    % AIRM(Ct, Ct_prev) — instantaneous covariance velocity

        % ===== Nonlinear eigenspace alignment =====
        subspaceAngleMax    % Max principal angle between Vt and V0 (degrees)
        subspaceAngleMean   % Mean principal angle between Vt and V0 (degrees)

        % ===== Threshold margin =====
        thresholdMarginMin  % min(threshold_dir - lambda) — negative = artifact
        thresholdMarginMean % mean margin across components
        nViolations         % number of components exceeding threshold

    end

    properties (Access = private)
        Ct_prev             % Previous Karcher-averaged covariance for step drift
    end

    methods

        function obj = riemannASRprobe(capacity)
            % RIEMANNASRPROBE - Construct probe with buffer capacity.
            %   capacity: max number of windows to store (default 1e5)
            if nargin < 1, capacity = 1e5; end
            obj = obj@asrProbe(capacity);
        end

        % ================================================================
        % setBaseline — override to read M directly (riemannASR has M, not vg0)
        % ================================================================
        function setBaseline(obj, asrObj)
            % SETBASELINE - Capture calibration state from a riemannASR instance.
            if isempty(asrObj.M)
                error('riemannASRprobe:notCalibrated', ...
                    'riemannASR must be calibrated before setBaseline().');
            end

            C0 = asrObj.M * asrObj.M';
            C0 = 0.5 * (C0 + C0');     % enforce symmetry

            [V0, D0]     = eig(C0, 'vector');
            [D0, idx]    = sort(D0, 'descend');
            V0           = V0(:, idx);

            obj.C0          = C0;
            obj.lambda0     = D0;
            obj.V0          = V0(:, 1:obj.subspaceRank);
            obj.T0          = asrObj.T;
            obj.C0_inv_sqrt = real(inv(sqrtm(C0)));
            obj.Ct_prev     = [];
        end

        % ================================================================
        % update — core per-window telemetry
        % ================================================================
        function update(obj, Ct, SCM, lambda, threshold_dir, R, trivial, riem_dist)
            % UPDATE - Record one window of telemetry.
            %   Ct            : Karcher-averaged covariance [C x C]
            %   SCM           : Raw sample covariance [C x C]
            %   lambda        : Sorted eigenvalues of Ct [C x 1]
            %   threshold_dir : Per-component threshold vector [C x 1]
            %   R             : Reconstruction matrix ([] if trivial)
            %   trivial       : logical
            %   riem_dist     : scalar AIRM distance from baseline

            if ~obj.enabled, return; end

            % Call base class update with the standard signature
            update@asrProbe(obj, lambda, threshold_dir, R, trivial, riem_dist);

            % ---- Karcher vs SCM distance ----
            obj.karcherVsScmDist(obj.step) = ...
                obj.calcAIRM(Ct, SCM);

            % ---- Step-to-step drift ----
            if ~isempty(obj.Ct_prev)
                obj.riemannDriftStep(obj.step) = ...
                    obj.calcAIRM(Ct, obj.Ct_prev);
            else
                obj.riemannDriftStep(obj.step) = 0;
            end
            obj.Ct_prev = Ct;

            % ---- Subspace principal angles (Vt vs V0) ----
            if ~isempty(obj.C0_inv_sqrt) && ~isempty(obj.V0)
                [Vt, ~] = eig(0.5*(Ct+Ct'), 'vector');
                Vt      = Vt(:, end:-1:end-obj.subspaceRank+1);  % top-k descending
                angles  = subspace(Vt, obj.V0) * (180/pi);
                obj.subspaceAngleMax(obj.step)  = max(angles);
                obj.subspaceAngleMean(obj.step) = mean(angles);
            else
                obj.subspaceAngleMax(obj.step)  = NaN;
                obj.subspaceAngleMean(obj.step) = NaN;
            end

            % ---- Threshold margin ----
            if ~isempty(threshold_dir) && ~isempty(lambda)
                margin = threshold_dir(:) - lambda(:);
                obj.thresholdMarginMin(obj.step)  = min(margin);
                obj.thresholdMarginMean(obj.step) = mean(margin);
                obj.nViolations(obj.step)         = sum(margin < 0);
            else
                obj.thresholdMarginMin(obj.step)  = NaN;
                obj.thresholdMarginMean(obj.step) = NaN;
                obj.nViolations(obj.step)         = 0;
            end
        end

        % ================================================================
        % reset — clear subclass buffers then call base reset
        % ================================================================
        function reset(obj)
            n = obj.subspaceRank;
            obj.karcherVsScmDist    = zeros(1, n);
            obj.riemannDriftStep    = zeros(1, n);
            obj.subspaceAngleMax    = zeros(1, n);
            obj.subspaceAngleMean   = zeros(1, n);
            obj.thresholdMarginMin  = zeros(1, n);
            obj.thresholdMarginMean = zeros(1, n);
            obj.nViolations         = zeros(1, n);
            obj.Ct_prev             = [];
            reset@asrProbe(obj);
        end

    end

    methods (Access = private)

        function d = calcAIRM(~, A, B)
            % CALCAIRM - Affine Invariant Riemannian Metric between two SPD matrices.
            %   d = || log(A^{-1/2} B A^{-1/2}) ||_F
            try
                A     = 0.5*(A + A');
                B     = 0.5*(B + B');
                As    = real(sqrtm(A));
                As_inv = inv(As);
                M_mid  = As_inv * B * As_inv;
                M_mid  = 0.5*(M_mid + M_mid');
                d      = norm(real(logm(M_mid)), 'fro');
            catch
                d = NaN;
            end
        end

    end
end
