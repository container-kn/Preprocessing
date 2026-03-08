classdef graphASRprobe < asrProbe
    % graphASRprobe - Diagnostic probe for graphASR.
    %   Extends asrProbe with:
    %     - Graph-spectral metrics (centroid, low/high energy ratio)
    %     - Participation ratio (effective dimensionality)
    %     - Step-to-step Riemannian drift
    %     - Cross-basis disagreement (optional)
    %
    %   Usage:
    %       p = graphASRprobe(12);
    %       p.setBaseline(graphASRobj);      % reads .vg0 and .lam
    %       p.update(lambda, thr, R, trivial, riem_dist);
    %       p.update(lambda, thr, R, trivial, riem_dist, specDis, subAng, commNorm);

    properties (Access = public)
        % ===== Graph Spectral Baseline =====
        lamGraph            % Graph Laplacian eigenvalues (from calibration)

        % ===== Spectral Shape =====
        participationRatio  % (sum lambda)^2 / sum(lambda^2) — effective dim
        spectralCentroid    % Weighted mean graph frequency
        lowHighEnergyRatio  % Energy in low vs high graph frequencies

        % ===== Step-to-Step Drift =====
        riemannDriftStep    % AIRM distance between t and t-1

        % ===== Cross-Basis Disagreement (Optional) =====
        spectralDisagreement
        subspaceAngle
        commutatorNorm

        % ===== Internal =====
        prevLambda          % Lambda from previous window for step drift
    end

    methods

        function obj = graphASRprobe(subspaceRank)
            if nargin < 1, subspaceRank = []; end
            obj = obj@asrProbe(subspaceRank);
            obj.resetOwn();
        end

        function setBaseline(obj, baselineObj)
            % SETBASELINE - Read graph-specific baseline from a graphASR object.
            %   Reads .vg0 (baseline per-frequency variance) and .lam (graph freqs).
            if ~obj.enabled, return; end

            % Build a minimal proxy so base setBaseline can work if M exists
            % For graphASR, there is no M — baseline is vg0 and lam directly.
            if isprop(baselineObj, 'vg0') && ~isempty(baselineObj.vg0)
                lambda0 = baselineObj.vg0(:);
                lambda0(lambda0 < 1e-12) = 1e-12;
                obj.lambda0      = sort(lambda0, 'descend');
                obj.subspaceRank = length(obj.lambda0);
            elseif isprop(baselineObj, 'M') && ~isempty(baselineObj.M)
                % Fallback: standard covariance baseline
                setBaseline@asrProbe(obj, baselineObj);
            else
                error('graphASRprobe:setBaseline', ...
                    'baselineObj must have property vg0 or M.');
            end

            if isprop(baselineObj, 'lam') && ~isempty(baselineObj.lam)
                obj.lamGraph = baselineObj.lam(:);
            end
        end

        function update(obj, lambda, threshold_dir, R, trivial, riem_dist, ...
                        specDis, subspaceAng, commNorm)
            % UPDATE - Record per-window metrics for graphASR.
            %   lambda        : [C x 1] per-frequency variance (graph domain)
            %   threshold_dir : [C x 1] per-component thresholds (or [])
            %   R             : [C x C] reconstruction matrix (or [])
            %   trivial       : logical scalar
            %   riem_dist     : scalar log-Euclidean distance from baseline
            %   specDis       : (optional) spectral disagreement scalar
            %   subspaceAng   : (optional) subspace angle scalar
            %   commNorm      : (optional) commutator norm scalar

            if ~obj.enabled, return; end

            if nargin < 7,  specDis      = NaN; end
            if nargin < 8,  subspaceAng  = NaN; end
            if nargin < 9,  commNorm     = NaN; end

            % 1. Base class shared metrics
            update@asrProbe(obj, lambda, threshold_dir, R, trivial, riem_dist);

            % Retrieve sanitised lambda from base (already sorted descending)
            lambda = obj.sanitiseLambda(lambda);

            % 2. Participation ratio
            obj.participationRatio(end+1,1) = sum(lambda)^2 / sum(lambda.^2);

            % 3. Graph spectral centroid and low/high energy ratio
            if ~isempty(obj.lamGraph)
                n       = min(length(obj.lamGraph), length(lambda));
                lg      = obj.lamGraph(1:n);
                lam_n   = lambda(1:n);
                obj.spectralCentroid(end+1,1) = ...
                    sum(lg .* lam_n) / sum(lam_n);
                mid = floor(n / 2);
                if mid > 0 && mid < n
                    obj.lowHighEnergyRatio(end+1,1) = ...
                        sum(lam_n(1:mid)) / sum(lam_n(mid+1:end));
                else
                    obj.lowHighEnergyRatio(end+1,1) = NaN;
                end
            else
                obj.spectralCentroid(end+1,1)    = NaN;
                obj.lowHighEnergyRatio(end+1,1)  = NaN;
            end

            % 4. Step-to-step Riemannian drift
            if ~isempty(obj.prevLambda)
                k = min(length(obj.prevLambda), length(lambda));
                obj.riemannDriftStep(end+1,1) = ...
                    norm(log(lambda(1:k) ./ obj.prevLambda(1:k)));
            else
                obj.riemannDriftStep(end+1,1) = NaN;
            end
            obj.prevLambda = lambda;

            % 5. Cross-basis disagreement (optional)
            obj.spectralDisagreement(end+1,1) = specDis;
            obj.subspaceAngle(end+1,1)        = subspaceAng;
            obj.commutatorNorm(end+1,1)       = commNorm;
        end

        function reset(obj)
            reset@asrProbe(obj);
            obj.resetOwn();
        end

    end

    methods (Access = private)
        function resetOwn(obj)
            obj.prevLambda           = [];
            obj.participationRatio   = zeros(0,1);
            obj.spectralCentroid     = zeros(0,1);
            obj.lowHighEnergyRatio   = zeros(0,1);
            obj.riemannDriftStep     = zeros(0,1);
            obj.spectralDisagreement = zeros(0,1);
            obj.subspaceAngle        = zeros(0,1);
            obj.commutatorNorm       = zeros(0,1);
        end
    end

end
