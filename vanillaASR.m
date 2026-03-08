classdef vanillaASR < handle
    % vanillaASR - Standard Artifact Subspace Reconstruction (ASR) Implementation.
    %   This class implements the classic ASR algorithm for real-time or offline
    %   EEG artifact removal. It identifies high-variance signal components
    %   relative to a calibrated baseline and reconstructs them using
    %   spherical subspace projection.
    %
    %   Reference:
    %       Mullen, T. R., et al. (2015). "Real-time neuroimaging and
    %       cognitive monitoring using wearable EEG." IEEE TBME.

    properties (Access = public)
        % Core Dimensions
        nchans              % Number of channels in the input signal
        nsamples = 0        % Total number of samples processed
        modifiedMask        % Boolean mask of samples where reconstruction occurred
        epsi = 1e-12        % FIX #4: was 10e-12 (= 1e-11), corrected to 1e-12

        % Algorithmic Parameters
        srate = 500         % Sampling rate in Hz
        cutoff = 20         % Standard deviation cutoff (k) for artifact detection
        blocksize = 10      % Block size for covariance estimation
        window_len = 0.5    % Window length (seconds) for moving statistics
        window_overlap = 0.66 % Overlap fraction between successive windows
        max_dropout_fraction = 0.1  % Parameter for EEG distribution fitting
        min_clean_fraction = 0.25   % Parameter for EEG distribution fitting
        maxdims = 0.66      % Max dimensions to reconstruct (fraction of nchans)
        lookahead = 0.25    % Lookahead buffer length in seconds
        stepsize = 32       % Samples between successive filter updates
        blockCount = 0      % FIX #12: renamed from 'step' to avoid MATLAB built-in conflict

        % Calibration State
        M                   % Mixing matrix (sqrt of baseline covariance)
        T                   % Threshold operator matrix
        V                   % Eigenvectors of the baseline covariance
        threshold_dir       % Per-component variance thresholds
        C0_inv_sqrt         % Inverse square root of baseline covariance (Riemannian)

        % Filter State (IIR Spectral Shaping)
        B                   % Numerator coefficients (Yule-Walker)
        A                   % Denominator coefficients (Yule-Walker)
        iir_state = []      % State of the spectral filter

        % Runtime & Buffering
        carry = []          % Boundary buffer to handle lookahead/causality
        last_R = []         % Reconstruction matrix from the previous step
        last_trivial = true % Boolean flag for previous step's reconstruction status

        % Diagnostic Probes
        probeRaw            % vanillaASRprobe instance for original signal statistics
        probeClean          % vanillaASRprobe instance for post-processed statistics
        timeProbe = []      % Object for high-resolution execution timing
    end

    methods

        function obj = vanillaASR(X, srate, cutoff_val, blocksize_val, timeProbe)
            % Construct an ASR object and initialize spectral filters.
            %   X            : Initial data matrix [chans x samples]
            %   srate        : Sampling frequency (Hz)
            %   cutoff_val   : SD threshold (default 20)
            %   blocksize_val: Covariance block size (default 10)
            %   timeProbe    : Timing probe object (always required by harness)

            % FIX #5: Input validation
            if nargin < 2 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('vanillaASR: srate must be a positive numeric scalar.');
            end
            if isempty(X) || ~isnumeric(X)
                error('vanillaASR: X must be a non-empty numeric matrix.');
            end

            if nargin >= 3 && ~isempty(cutoff_val),   obj.cutoff    = cutoff_val;    end
            if nargin >= 4 && ~isempty(blocksize_val), obj.blocksize = blocksize_val; end
            if nargin >= 5,                             obj.timeProbe = timeProbe;     end

            obj.srate = srate;
            X = obj.determineSignalShape(X);
            [obj.nchans, obj.nsamples] = size(X);
            obj.modifiedMask = false(1, obj.nsamples);
            [obj.B, obj.A] = vanillaASR.getSpectralFilter(srate);
        end

        function X = determineSignalShape(obj, X)
            % FIX #6: After calibration, validate orientation against known nchans.
            %   Before calibration (nchans empty), fall back to the heuristic.
            if ~isnumeric(X) || isempty(X)
                error('vanillaASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            if ~isempty(obj.nchans)
                % Post-calibration: trust nchans to enforce correct orientation
                if size(X, 1) == obj.nchans
                    return;
                elseif size(X, 2) == obj.nchans
                    X = X';
                else
                    error('vanillaASR:channelMismatch', ...
                        'Data has %d rows and %d cols; expected %d channels.', ...
                        size(X,1), size(X,2), obj.nchans);
                end
            else
                % Pre-calibration: heuristic — more columns than rows = [chans x samples]
                if size(X, 1) > size(X, 2)
                    X = X';
                end
            end
        end

        function calibrate(obj, X)
            % Calibrate the ASR model using clean baseline data.
            % Calculates the baseline covariance structure and sets thresholds.

            % FIX #5: Validate input before proceeding
            if isempty(X) || ~isnumeric(X)
                error('vanillaASR:calibrate', 'Calibration data must be a non-empty numeric matrix.');
            end

            X = obj.determineSignalShape(X);

            % 1. Spectral filtering (shaping the baseline)
            obj.timeProbe.start('prefilter');
            X = filter(obj.B, obj.A, double(X), [], 2);
            obj.timeProbe.stop('prefilter');

            % 2. Robust Covariance Estimation
            % FIX #7: blockCovariance normalization corrected (see static method)
            U = vanillaASR.blockCovariance(X, obj.blocksize);
            U = (U + U') / 2;                          % Enforce symmetry
            U = U + obj.epsi * eye(obj.nchans);        % Tikhonov regularisation

            % 3. Subspace Decomposition
            obj.timeProbe.start('subspaceM');
            obj.M = sqrtm(real(U));
            obj.timeProbe.stop('subspaceM');

            % 4. Riemannian Metric Precomputation
            % Precompute C0^{-1/2} for Affine-Invariant Riemannian Distance (AIRM)
            [V0, D0] = eig(obj.M * obj.M');
            lambda0  = diag(D0);
            lambda0(lambda0 < obj.epsi) = obj.epsi;
            obj.C0_inv_sqrt = V0 * diag(1 ./ sqrt(lambda0)) * V0';

            % 5. Eigenvector basis for threshold projection
            obj.timeProbe.start('eigCal');
            [obj.V, ~] = eig(obj.M);
            obj.timeProbe.stop('eigCal');

            % 6. Threshold Estimation
            [obj.T, obj.threshold_dir] = obj.getThresholdVector(X);

            % 7. Initialise Diagnostic Probes
            obj.probeRaw   = vanillaASRprobe(12);
            obj.probeClean = vanillaASRprobe(12);
            obj.probeRaw.setBaseline(obj);
            obj.probeClean.setBaseline(obj);

            % Reset runtime state so process() starts clean
            obj.carry        = [];
            obj.last_R       = [];
            obj.last_trivial = true;
            obj.iir_state    = [];
            obj.blockCount   = 0;
        end

        function cleanedSignal = process(obj, data)
            % process - Reconstructs EEG signal via subspace projection with telemetry.

            if isempty(obj.M)
                error('vanillaASR:notCalibrated', 'Run calibrate() before process().');
            end

            % FIX #5: validate input dimensions
            if isempty(data) || ~isnumeric(data)
                error('vanillaASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            data = obj.determineSignalShape(data);  % FIX #6: uses nchans-aware version

            [C, S] = size(data);
            P = round(obj.lookahead * obj.srate);

            % --- Boundary & Buffer Management ---
            if isempty(obj.carry)
                % Mirror-reflect the first P samples for causal boundary handling
                mirror_idx   = 1 + mod(((P+1):-1:2) - 1, S);
                obj.carry    = repmat(2*data(:,1), 1, P) - data(:, mirror_idx);
            end

            data = [obj.carry, double(data)];
            data(~isfinite(data)) = 0;

            % FIX #8: Extend modifiedMask to cover the new samples being processed
            obj.modifiedMask(obj.nsamples + (1:S)) = false;

            % --- Spectral Shaping (Yule-Walker IIR) ---
            obj.timeProbe.start('spectral_shaping');
            [X, obj.iir_state] = filter(obj.B, obj.A, data, obj.iir_state, 2);
            obj.timeProbe.stop('spectral_shaping');

            N        = round(obj.window_len * obj.srate);
            update_at = obj.stepsize : obj.stepsize : size(X, 2);

            if isempty(obj.last_R)
                obj.last_R = eye(C);
                update_at  = [1, update_at];
            end

            last_n = 0;

            for j = 1:length(update_at)
                obj.blockCount = obj.blockCount + 1;
                n   = update_at(j);
                idx = max(1, n-N+1) : n;
                Xw  = X(:, idx);

                % 1. Local Covariance Estimation
                obj.timeProbe.start('local_covariance');
                Ct = (Xw * Xw') / size(Xw, 2);        % zero-mean covariance proxy
                Ct = (Ct + Ct') / 2 + 1e-10 * eye(C);
                obj.timeProbe.stop('local_covariance');

                % 2. Subspace Analysis — single sort, consistent ascending order
                % FIX #3: derive both lambda_asc and Vw_asc from ONE joint sort
                obj.timeProbe.start('subspace_eig');
                [Vw, Dmat]          = eig(Ct);
                lambda_raw          = diag(Dmat);
                lambda_raw(lambda_raw <= 0) = eps;
                [lambda_asc, asc_ord] = sort(lambda_raw, 'ascend');
                Vw_asc              = Vw(:, asc_ord);   % eigenvectors in matching order
                [lambda_desc, desc_ord] = sort(lambda_raw, 'descend');
                Vw_desc             = Vw(:, desc_ord);  % kept for probe reporting
                obj.timeProbe.stop('subspace_eig');

                % 3. Affine-Invariant Riemannian Distance (AIRM)
                obj.timeProbe.start('riemannian_distance');
                delta     = obj.C0_inv_sqrt * Ct * obj.C0_inv_sqrt;
                delta     = (delta + delta') / 2;
                lam_delta = eig(delta);
                riem_dist = sqrt(sum(log(max(lam_delta, 1e-12)).^2));
                obj.timeProbe.stop('riemannian_distance');

                % 4. Threshold Evaluation & Component Selection
                threshold_dir_local = sum((obj.T * Vw_asc).^2);

                keep    = (lambda_asc < threshold_dir_local') | ...
                          ((1:C)' <= C - round(obj.maxdims * C));
                trivial = all(keep);

                % 5. Reconstruction Matrix
                obj.timeProbe.start('reconstruction_matrix');
                if ~trivial
                    R = real(obj.M * pinv(bsxfun(@times, keep', Vw_asc' * obj.M)) * Vw_asc');
                else
                    R = eye(C);
                end
                obj.timeProbe.stop('reconstruction_matrix', ...
                    struct('trivial', trivial, 'dims_rejected', sum(~keep)));

                % 6. Signal Blending & Probe Updates
                if ~trivial || ~obj.last_trivial
                    subrange     = (last_n+1) : n;
                    raw_segment  = data(:, subrange);

                    % FIX #2: offset mask by P to account for the carry-prepend
                    mask_idx = obj.nsamples + subrange - P;
                    mask_idx = mask_idx(mask_idx >= 1);
                    obj.modifiedMask(mask_idx) = true;

                    blend = (1 - cos(pi * (1:length(subrange)) / length(subrange))) / 2;
                    cleaned_segment = bsxfun(@times, blend,   R           * raw_segment) + ...
                                      bsxfun(@times, 1-blend, obj.last_R  * raw_segment);
                    data(:, subrange) = cleaned_segment;

                    % FIX #1 & #9: update both probes; probeClean uses the cleaned segment
                    obj.probeRaw.update(Ct, Vw_desc, lambda_desc, threshold_dir_local, ...
                                        R, trivial, riem_dist, obj.T);

                    Ct_clean   = (cleaned_segment * cleaned_segment') / size(cleaned_segment, 2);
                    Ct_clean   = (Ct_clean + Ct_clean') / 2 + 1e-10 * eye(C);
                    [Vc, Dc]   = eig(Ct_clean);
                    lam_c      = diag(Dc);
                    [lam_c_d, cord] = sort(lam_c, 'descend');
                    obj.probeClean.update(Ct_clean, Vc(:,cord), lam_c_d, threshold_dir_local, ...
                                          R, trivial, riem_dist, obj.T);
                end

                obj.last_R       = R;
                obj.last_trivial = trivial;
                last_n           = n;
            end

            % --- Output Management ---
            obj.carry        = data(:, end-P+1 : end);
            cleanedSignal    = data(:, P+1 : end-P);   % strip carry prefix AND lookahead suffix
            obj.nsamples     = obj.nsamples + S;
        end

    end % end public methods


    methods (Static)

        function [B, A] = getSpectralFilter(srate)
            % Design a Yule-Walker IIR filter for EEG spectral shaping.
            nyq      = srate / 2;
            freqvals = [0 2 3 13 16 40 min(80, nyq-1)] * 2 / srate;
            freqvals = [freqvals, 1];
            amps     = [3 0.75 0.33 0.33 1 1 3 3];
            [B, A]   = yulewalk(8, freqvals, amps);
        end

        function U = blockCovariance(X, blocksize)
            % Compute a block-averaged covariance estimate.
            %
            % FIX #7: Original code divided by length(1:blocksize:S) — the number
            % of starting offsets — while accumulating over k = 1:blocksize offsets.
            % These counts differ, producing a biased estimate. The fix accumulates
            % all sample outer-products and divides by the actual sample count.
            %
            % Method: stride through X in blocks of `blocksize`, accumulate X*X',
            % then normalise by total number of columns used.
            [~, S]    = size(X);
            U         = zeros(size(X, 1));
            n_used    = 0;

            for k = 1:blocksize
                idx = k : blocksize : S;
                if isempty(idx), continue; end
                Xk     = X(:, idx);
                U      = U + (Xk * Xk');
                n_used = n_used + size(Xk, 2);
            end

            if n_used == 0
                error('vanillaASR:blockCovariance', 'No samples available for covariance estimation.');
            end
            U = U / n_used;
        end

    end % end static methods


    methods

        function [T, threshold_dir] = getThresholdVector(obj, X)
            % Estimate the per-component threshold from the baseline distribution.
            % Uses a truncated Gaussian fit (fit_eeg_distribution) on windowed RMS.
            [C, S] = size(X);
            N      = round(obj.window_len * obj.srate);
            Xp     = abs(obj.V' * X);   % Project to eigenvector space

            mu  = zeros(1, C);
            sig = zeros(1, C);

            for c = 1:C
                rms_sq = Xp(c, :).^2;

                % Build matrix of windowed RMS values
                starts  = 1 : round(N * (1 - obj.window_overlap)) : S - N;
                if isempty(starts)
                    warning('vanillaASR:getThresholdVector', ...
                        'Channel %d: too few samples for windowed RMS; using full signal.', c);
                    mu(c)  = sqrt(mean(rms_sq));
                    sig(c) = std(sqrt(rms_sq));
                    continue;
                end
                idx_mat     = bsxfun(@plus, starts(:), (0:N-1));
                idx_mat     = min(idx_mat, S);               % clamp to signal length
                rms_windowed = sqrt(sum(rms_sq(idx_mat), 2)' / N);

                [mu(c), sig(c)] = fit_eeg_distribution(rms_windowed, ...
                    obj.min_clean_fraction, obj.max_dropout_fraction);
            end

            threshold_vals = mu + obj.cutoff * sig;
            T              = diag(threshold_vals) * obj.V';
            threshold_dir  = threshold_vals .^ 2;
        end

        function stepCount = getCurrentStep(obj)
            % Return current block iteration count.
            % FIX #12: property renamed from step -> blockCount
            stepCount = obj.blockCount;
        end

    end

end
