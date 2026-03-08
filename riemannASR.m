classdef riemannASR < handle
    % riemannASR - Riemannian Geometry Artifact Subspace Reconstruction.
    %   Extends vanilla ASR by replacing the per-window sample covariance
    %   with a geodesically-smoothed estimate computed as the Karcher mean
    %   of the current SCM and the previous window's covariance on the SPD
    %   manifold. This makes artifact detection more robust to transient
    %   covariance spikes without requiring an explicit adaptation model.
    %
    %   Key differences from vanillaASR:
    %     - Covariance smoothed via positive_definite_karcher_mean() each window
    %     - Eigenspace computed by rasr_nonlinear_eigenspace() (manifold-aware)
    %     - Spectral shaping via Yule-Walker IIR filter (same as RASR)
    %
    %   External dependencies (must be on MATLAB path):
    %     block_geometric_median  — RASR toolbox
    %     fit_eeg_distribution    — EEGLAB / RASR toolbox
    %   (rasr_nonlinear_eigenspace and positive_definite_karcher_mean are no
    %    longer required — replaced with regularised eig and inline SPD midpoint)
    %
    %   Constructor signature (matches Experimenting.m harness):
    %       obj = riemannASR(srate, params)
    %       obj = riemannASR(srate, params, tProbe)

    properties (Access = public)

        %% Core
        nchans
        srate
        nsamples    = 0
        modifiedMask

        %% Parameters
        cutoff               = 5       % SD threshold multiplier
        blocksize            = 10      % Covariance block size (samples) for calibration
        window_len           = 0.5     % Window length (seconds) for threshold estimation
        window_overlap       = 0.66    % Overlap fraction for windowed RMS
        max_dropout_fraction = 0.1     % fit_eeg_distribution parameter
        min_clean_fraction   = 0.25    % fit_eeg_distribution parameter
        maxdims              = 0.66    % Max fraction of dims to reconstruct per window
        lookahead            = 0.25    % Lookahead buffer length (seconds)
        stepsize             = 32      % Reconstruction step size (samples)

        %% Calibration state
        M               % Mixing matrix: sqrtm of geometric median covariance
        T               % Threshold matrix: diag(mu + cutoff*sig) * V'
        iir_B           % IIR filter numerator (spectral shaping)
        iir_A           % IIR filter denominator

        %% Runtime state
        carry           = []       % Lookahead carry buffer
        iir_state       = []       % IIR filter state (persists across process calls)
        last_R          = []       % Last reconstruction matrix
        last_trivial    = true     % Trivial flag from last window
        cov_prev        = []       % Previous Karcher-averaged covariance

        %% Telemetry
        timeProbe                  % timeProbe instance for stage timing
        probeRaw                   % riemannASRprobe on raw signal
        probeClean                 % riemannASRprobe on cleaned signal

    end

    methods

        %% ================= CONSTRUCTOR =================
        function obj = riemannASR(srate, params, tProbe)
            % RIEMANNASR - Construct a riemannASR instance.
            %   srate  : Sampling rate in Hz (required)
            %   params : Struct of parameter overrides (optional)
            %   tProbe : timeProbe instance (optional — created internally if omitted)

            if nargin < 1 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('riemannASR:invalidSrate', ...
                    'srate must be a positive numeric scalar.');
            end

            obj.srate = srate;

            % Build spectral shaping filter
            [obj.iir_B, obj.iir_A] = riemannASR.buildSpectralFilter(srate);

            % timeProbe
            if nargin >= 3 && ~isempty(tProbe)
                obj.timeProbe = tProbe;
            else
                obj.timeProbe = timeProbe();
            end

            % Apply parameter overrides
            if nargin >= 2 && isstruct(params)
                fields = fieldnames(params);
                for i = 1:numel(fields)
                    if isprop(obj, fields{i})
                        obj.(fields{i}) = params.(fields{i});
                    end
                end
            end

            obj.modifiedMask = false(1, 0);
        end


        %% ================= CALIBRATION =================
        function calibrate(obj, X)
            % CALIBRATE - Estimate mixing matrix M and threshold T from clean data.
            %   X : [nchans x ntime] calibration signal

            if isempty(X) || ~isnumeric(X)
                error('riemannASR:calibrate', ...
                    'Calibration data must be a non-empty numeric matrix.');
            end
            if size(X, 1) > size(X, 2)
                X = X';
            end

            fprintf('--- [riemannASR] Calibrating ---\n');
            obj.timeProbe.start('calibration');

            [obj.nchans, ~] = size(X);

            % Spectral shaping — initialise filter state from calibration data
            [X_filt, obj.iir_state] = filter(obj.iir_B, obj.iir_A, double(X), [], 2);
            X_filt = X_filt';   % [ntime x nchans] for block operations

            % ---- Block geometric median covariance ----
            nBlocks = length(1:obj.blocksize:size(X_filt, 1));
            U       = zeros(nBlocks, obj.nchans^2);
            idx     = 1:obj.blocksize:size(X_filt, 1);

            for k = 1:numel(idx)
                rng  = idx(k) : min(size(X_filt,1), idx(k)+obj.blocksize-1);
                blk  = X_filt(rng, :);
                U(k,:) = reshape(blk' * blk, 1, []);
            end

            Cmed   = reshape(block_geometric_median(U / obj.blocksize), ...
                             obj.nchans, obj.nchans);
            obj.M  = sqrtm(real(Cmed));

            % ---- Eigenspace for threshold estimation (regularised eig) ----
            Cmed_reg = riemannASR.regulariseSPD(Cmed);
            [V, Dmat] = eig(Cmed_reg);
            [~, ord]  = sort(real(diag(Dmat)));
            V         = real(V(:, ord));

            % ---- Windowed RMS threshold estimation ----
            Xproj = abs(X_filt * V);
            N     = round(obj.window_len * obj.srate);
            mu    = zeros(1, obj.nchans);
            sig   = zeros(1, obj.nchans);

            for c = 1:obj.nchans
                rms_sq  = Xproj(:, c) .^ 2;
                idx_win = round(1 : N*(1-obj.window_overlap) : size(Xproj,1)-N);
                idx_mat = bsxfun(@plus, idx_win', (0:N-1));
                rmsw    = sqrt(sum(rms_sq(idx_mat), 2) / N);

                [mu(c), sig(c)] = fit_eeg_distribution( ...
                    rmsw, obj.min_clean_fraction, obj.max_dropout_fraction);
            end

            obj.T = diag(mu + obj.cutoff * sig) * V';

            % ---- Initialise probes ----
            obj.probeRaw   = riemannASRprobe(12);
            obj.probeClean = riemannASRprobe(12);
            obj.probeRaw.setBaseline(obj);
            obj.probeClean.setBaseline(obj);

            % ---- Reset runtime state ----
            obj.nsamples     = 0;
            obj.modifiedMask = false(1, 0);
            obj.carry        = [];
            obj.last_R       = [];
            obj.last_trivial = true;
            obj.cov_prev     = [];

            obj.timeProbe.stop('calibration');
            fprintf('--- [riemannASR] Calibration done ---\n');
        end


        %% ================= PROCESS =================
        function cleaned = process(obj, data)
            % PROCESS - Clean artifact-contaminated EEG using Riemannian ASR.
            %   data : [nchans x ntime] input signal
            %   cleaned : [nchans x ntime] output signal (same size)

            if isempty(obj.M)
                error('riemannASR:notCalibrated', ...
                    'Call calibrate() before process().');
            end
            if size(data, 1) > size(data, 2)
                data = data';
            end
            if size(data, 1) ~= obj.nchans
                error('riemannASR:channelMismatch', ...
                    'Expected %d channels, got %d.', obj.nchans, size(data,1));
            end

            obj.timeProbe.start('process');

            [C, S] = size(data);
            P      = round(obj.lookahead * obj.srate);

            % ---- Lookahead carry buffer ----
            if isempty(obj.carry)
                obj.carry = repmat(2*data(:,1), 1, P) ...
                    - data(:, 1 + mod(((P+1):-1:2)-1, S));
            end

            data = [obj.carry, data];
            data(~isfinite(data)) = 0;

            % ---- Spectral shaping (continues from calibration iir_state) ----
            [X, obj.iir_state] = filter(obj.iir_B, obj.iir_A, ...
                                        double(data), obj.iir_state, 2);

            N         = round(obj.window_len * obj.srate);
            update_at = obj.stepsize : obj.stepsize : size(X, 2);

            if isempty(obj.last_R)
                update_at  = [1, update_at];
                obj.last_R = eye(C);
            end

            % Pre-extend modifiedMask for this block (exclude lookahead prefix)
            nOut = S;
            obj.modifiedMask(obj.nsamples + (1:nOut)) = false;

            last_n = 0;

            for j = 1:numel(update_at)

                n   = update_at(j);
                idx = max(1, n-N+1) : n;
                Xw  = X(:, idx);

                % ---- Sample covariance ----
                SCM = (1/numel(idx)) * (Xw * Xw');
                SCM = 0.5 * (SCM + SCM');   % enforce symmetry

                % ---- Riemannian geodesic midpoint with previous covariance ----
                % Replaces positive_definite_karcher_mean (RASR toolbox dependency).
                % For exactly 2 SPD matrices A, B the geodesic midpoint on the
                % SPD manifold is: A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}
                if ~isempty(obj.cov_prev)
                    Ct = riemannASR.spdGeodesicMidpoint(obj.cov_prev, SCM);
                else
                    Ct = SCM;
                end

                obj.cov_prev = Ct;

                % ---- Eigenspace (regularised eig, replaces rasr_nonlinear_eigenspace) ----
                % rasr_nonlinear_eigenspace frequently hits max iterations and
                % returns near-singular results on real EEG data. Standard eig
                % on a Tikhonov-regularised SPD covariance is stable and faster.
                Ct_reg   = riemannASR.regulariseSPD(Ct);
                [V, Dmat] = eig(Ct_reg);
                [D, ord] = sort(real(diag(Dmat)));
                V        = real(V(:, ord));

                % ---- Thresholding ----
                threshold_dir = sum((obj.T * V).^2)';
                keep    = D < threshold_dir | ...
                          (1:C)' <= C - round(obj.maxdims * C);
                trivial = all(keep);

                % ---- Reconstruction matrix ----
                if ~trivial
                    R = real(obj.M * ...
                        pinv(bsxfun(@times, keep', V' * obj.M)) * V');
                else
                    R = eye(C);
                end

                % ---- Overlap-add blending ----
                if ~trivial || ~obj.last_trivial
                    subrange = (last_n+1) : n;
                    blend    = (1 - cos(pi*(1:numel(subrange)) / numel(subrange))) / 2;

                    data(:, subrange) = ...
                        bsxfun(@times, blend,   R            * data(:, subrange)) + ...
                        bsxfun(@times, 1-blend, obj.last_R   * data(:, subrange));

                    % Mark reconstructed samples in modifiedMask
                    % Offset by P (lookahead prefix length) to align with output
                    out_range = subrange - P;
                    valid     = out_range >= 1 & out_range <= nOut;
                    obj.modifiedMask(obj.nsamples + out_range(valid)) = true;
                end

                % ---- Probe update ----
                riem_dist = obj.calcBaselineDist(Ct);

                obj.probeRaw.update( ...
                    Ct, SCM, D, threshold_dir, ...
                    R, trivial, riem_dist);

                if ~trivial || ~obj.last_trivial
                    % Clean probe: read from data AFTER blend
                    Xw_clean    = data(:, idx);
                    Xw_clean_f  = filter(obj.iir_B, obj.iir_A, ...
                                         double(Xw_clean), [], 2);
                    SCM_clean   = (1/numel(idx)) * (Xw_clean_f * Xw_clean_f');
                    SCM_clean   = 0.5*(SCM_clean + SCM_clean');
                    SCM_clean_reg = riemannASR.regulariseSPD(SCM_clean);
                    [~, Dc_mat]   = eig(SCM_clean_reg);
                    Dc            = sort(real(diag(Dc_mat)));
                    riem_clean  = obj.calcBaselineDist(SCM_clean);

                    obj.probeClean.update( ...
                        SCM_clean, SCM_clean, Dc, threshold_dir, ...
                        R, trivial, riem_clean);
                end

                obj.last_R       = R;
                obj.last_trivial = trivial;
                last_n           = n;
            end

            obj.carry    = data(:, end-P+1 : end);
            cleaned      = data(:, P+1 : end);    % strip lookahead prefix
            obj.nsamples = obj.nsamples + nOut;

            obj.timeProbe.stop('process');
        end


        %% ================= RESET =================
        function reset(obj)
            % RESET - Clear runtime state between independent processing runs.
            %   Does NOT reset calibration (M, T) or filter coefficients.
            %   Call this between unrelated processing segments if you don't
            %   want carry/iir_state contamination across them.
            obj.carry        = [];
            obj.iir_state    = [];
            obj.last_R       = [];
            obj.last_trivial = true;
            obj.cov_prev     = [];
            obj.nsamples     = 0;
            obj.modifiedMask = false(1, 0);
        end

    end % public methods


    methods (Access = private)

        function d = calcBaselineDist(obj, Ct)
            % CALCBASELINEDIST - AIRM distance from calibration covariance C0.
            %   Uses eigendecomposition-based log to avoid logm warnings on
            %   matrices with near-zero or negative eigenvalues.
            try
                C0  = riemannASR.regulariseSPD(obj.M * obj.M');
                Ctr = riemannASR.regulariseSPD(Ct);
                As     = real(sqrtm(C0));
                As_inv = pinv(As);
                M_mid  = As_inv * Ctr * As_inv;
                M_mid  = riemannASR.regulariseSPD(M_mid);
                % Eigendecomposition-based matrix log — avoids logm on
                % matrices with numerically non-positive eigenvalues
                [Vl, Dl] = eig(M_mid);
                dl       = real(diag(Dl));
                dl       = max(dl, 1e-12);
                logM     = real(Vl * diag(log(dl)) * Vl');
                d        = norm(logM, 'fro');
            catch
                d = NaN;
            end
        end

    end % private methods


    methods (Static, Access = private)

        function C = regulariseSPD(C)
            % REGULARISESPD - Tikhonov regularisation to ensure strict SPD.
            %   Adds epsilon * trace(C)/n * I to push all eigenvalues positive.
            C   = 0.5 * (C + C');
            n   = size(C, 1);
            eps = max(1e-6 * trace(C) / n, 1e-10);
            C   = C + eps * eye(n);
        end

        function Cm = spdGeodesicMidpoint(A, B)
            % SPDGEODESICMIDPOINT - Geodesic midpoint of two SPD matrices.
            %   Replaces positive_definite_karcher_mean for the 2-matrix case.
            %   Formula: A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}
            A  = riemannASR.regulariseSPD(A);
            B  = riemannASR.regulariseSPD(B);
            As = real(sqrtm(A));
            Ai = pinv(As);
            M  = Ai * B * Ai;
            M  = riemannASR.regulariseSPD(M);
            Cm = real(As * sqrtm(M) * As);
            Cm = 0.5 * (Cm + Cm');
        end

        function [B, A] = buildSpectralFilter(srate)
            % BUILDSPECTRALFILTER - Yule-Walker IIR spectral shaping filter.
            %   Same design as RASR / original ASR — emphasises artifact-band
            %   frequencies and de-emphasises low/high ends.
            nyq      = srate / 2;
            freqvals = [0 2 3 13 16 40 min(80, nyq-1)] * 2 / srate;
            freqvals = [freqvals 1];
            amps     = [3 0.75 0.33 0.33 1 1 3 3];
            [B, A]   = yulewalk(8, freqvals, amps);
        end

    end % static methods

end
