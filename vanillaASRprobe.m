classdef vanillaASRprobe < asrProbe
    % vanillaASRprobe - Diagnostic probe for vanillaASR and originalASR.
    %   Extends asrProbe with:
    %     - Subspace geometry (principal angle between current and baseline V)
    %     - Threshold drift (how much T has moved from calibration)
    %
    %   Usage:
    %       p = vanillaASRprobe(12);
    %       p.setBaseline(asrObj);           % asrObj must have .M and .T
    %       p.update(Ct, Vc, lambda, thr, R, trivial, riem_dist, T_current);

    properties (Access = public)
        % ===== Subspace Geometry =====
        subspaceAngleMax    % max principal angle between V_t and V_0
        subspaceAngleMean   % mean principal angle

        % ===== Threshold Drift =====
        thresholdDrift      % ||diag(T_t) - T_0|| / ||T_0||
    end

    methods

        function obj = vanillaASRprobe(subspaceRank)
            if nargin < 1, subspaceRank = []; end
            obj = obj@asrProbe(subspaceRank);   % call base constructor
            obj.resetOwn();
        end

        function update(obj, Ct, Vc, lambda, threshold_dir, ...
                        R, trivial, riem_dist, T_current)
            % UPDATE - Record all metrics for vanillaASR / originalASR.
            %   Ct            : [C x C] current window covariance
            %   Vc            : [C x C] current eigenvectors (descending order)
            %   lambda        : [C x 1] current eigenvalues (descending)
            %   threshold_dir : [C x 1] per-component thresholds
            %   R             : [C x C] reconstruction matrix
            %   trivial       : logical scalar
            %   riem_dist     : scalar AIRM distance
            %   T_current     : [C x C] current threshold operator T

            if ~obj.enabled, return; end

            % 1. Base class shared metrics
            update@asrProbe(obj, lambda, threshold_dir, R, trivial, riem_dist);

            % 2. Subspace angles (requires Vc and baseline V0)
            if ~isempty(obj.V0) && ~isempty(Vc)
                k      = min(obj.subspaceRank, size(Vc, 2));
                s      = svd(obj.V0' * Vc(:, 1:k));
                s      = min(max(s, -1), 1);
                angles = acos(s);
                obj.subspaceAngleMax(end+1,1)  = max(angles);
                obj.subspaceAngleMean(end+1,1) = mean(angles);
            else
                obj.subspaceAngleMax(end+1,1)  = NaN;
                obj.subspaceAngleMean(end+1,1) = NaN;
            end

            % 3. Threshold drift
            if ~isempty(obj.T0) && ~isempty(T_current)
                Tt = diag(T_current);
                obj.thresholdDrift(end+1,1) = norm(Tt - obj.T0) / norm(obj.T0);
            else
                obj.thresholdDrift(end+1,1) = NaN;
            end

            % 4. Energy ratio using full Ct (more accurate than lambda-only version)
            if ~isempty(obj.V0) && ~isempty(Ct)
                % Override base class energyRatio with the full-matrix version
                obj.energyRatio(end,1) = trace(obj.V0' * Ct * obj.V0) / trace(Ct);
            end
        end

        function reset(obj)
            reset@asrProbe(obj);    % clear base metrics
            obj.resetOwn();
        end

    end

    methods (Access = private)
        function resetOwn(obj)
            obj.subspaceAngleMax  = zeros(0,1);
            obj.subspaceAngleMean = zeros(0,1);
            obj.thresholdDrift    = zeros(0,1);
        end
    end

end
