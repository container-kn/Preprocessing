classdef emaASR < handle
    % emaASR - Adaptive Artifact Subspace Reconstruction with EMA & Telemetry.
    %   This implementation utilizes Exponential Moving Averages (EMA) to adapt
    %   baseline statistics in real-time based on a rolling buffer of clean EEG.
    %
    %   Constructor signature (matches Experimenting.m harness):
    %       obj = emaASR(srate, params)
    %       obj = emaASR(srate, params, tProbe)

    properties (Access = public)
        % --- Core Dimensions ---
        nchans
        nsamples = 0
        srate = 500

        % --- ASR Parameters ---
        cutoff = 20
        window_len = 0.5
        lookahead = 0.25
        stepsize = 32
        maxdims = 0.66
        blocksize = 10

        % --- EMA Adaptation Parameters ---
        beta1 = 0.02          % Update rate for Mixing Matrix (M)
        beta2 = 0.03          % Update rate for Threshold Operator (T)
        cleanBuffLength = 150 % Rolling buffer length in seconds

        % --- State Variables ---
        state                 % Internal ASR state (carry, iir, cov, etc.)
        M                     % Current EMA-adapted Mixing Matrix
        T                     % Current EMA-adapted Threshold Operator

        % --- Buffers & Diagnostics ---
        cleanBuffer = []      % Rolling window of clean segments (starts empty)
        modifiedMask          % Boolean mask of reconstructed samples
        tProbe                % High-resolution timing utility
        emaProbe              % Diagnostic telemetry for adaptation
    end

    methods

        function obj = emaASR(srate, params, tProbe)
            % Constructor
            %   srate   : Sampling rate in Hz
            %   params  : Struct of algorithm parameters (fields mapped to properties)
            %   tProbe  : (optional) timeProbe instance — created if not supplied

            if nargin < 1 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('emaASR: srate must be a positive numeric scalar.');
            end
            obj.srate = srate;

            % FIX #10: fallback uses timeProbe (the actual class), not asrTimeProbe
            if nargin >= 3 && ~isempty(tProbe)
                obj.tProbe = tProbe;
            else
                obj.tProbe = timeProbe();
            end

            % Populate properties from the params struct
            if nargin >= 2 && isstruct(params)
                f = fieldnames(params);
                for i = 1:length(f)
                    if isprop(obj, f{i})
                        obj.(f{i}) = params.(f{i});
                    end
                end
            end

            obj.modifiedMask = false(1, 0);
        end

        function calibrate(obj, X)
            % CALIBRATE - Establish initial baseline with strict type enforcement.
            fprintf('--- [emaASR] Initializing Adaptive Baseline ---\n');

            if isempty(X) || ~isnumeric(X)
                error('emaASR:calibrate', 'Calibration data must be a non-empty numeric matrix.');
            end

            [obj.nchans, ~] = size(X);
            obj.emaProbe    = emaASRprobe(obj.nchans);

            fs_val     = double(obj.srate);
            cutoff_val = double(obj.cutoff);
            % FIX #13: win_len_val was computed but never used — removed.
            %          obj.window_len is passed directly below.

            obj.tProbe.start('calibration');
            obj.state = emaASR.calibrateBaseline( ...
                double(X), fs_val, cutoff_val, obj.blocksize, obj.window_len);
            obj.tProbe.stop('calibration');

            % Snapshot calibration M/T into obj for EMA tracking
            obj.M = obj.state.M;
            obj.T = obj.state.T;

            % Anchor the probe to the frozen calibration state
            obj.emaProbe.setBaseline(obj);

            % Initialise runtime state fields
            % NOTE: obj.state.iir is intentionally NOT reset here.
            % calibrateBaseline() saves the final IIR filter state from the
            % calibration pass into state.iir. Preserving it avoids a startup
            % transient on the first process() call.
            obj.state.carry        = [];
            obj.state.last_R       = [];
            obj.state.last_trivial = true;
            obj.state.cov          = [];

            % FIX #9: cleanBuffer starts empty — NOT pre-filled with zeros.
            % Pre-filling with zeros biases M_new toward zero until the buffer
            % fills. The rolling-append logic in adaptBaseline handles growth.
            obj.cleanBuffer = [];

            % Reset sample counter on re-calibration
            obj.nsamples     = 0;
            obj.modifiedMask = false(1, 0);
        end

        function cleanedSignal = process(obj, data)
            % PROCESS - Subspace reconstruction with adaptive baseline tracking.

            if isempty(obj.M)
                error('emaASR:notCalibrated', 'Run calibrate() before process().');
            end
            if isempty(data) || ~isnumeric(data)
                error('emaASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end

            [C, S] = size(data);
            P = round(obj.lookahead * obj.srate);
            N = round(obj.window_len * obj.srate);

            % 1. Boundary Handling
            if isempty(obj.state.carry)
                mirror_idx      = 1 + mod(((P+1):-1:2) - 1, S);
                obj.state.carry = repmat(2*data(:,1), 1, P) - data(:, mirror_idx);
            end
            data = [obj.state.carry, double(data)];
            data(~isfinite(data)) = 0;

            % FIX #1 (partial): extend modifiedMask before the loop
            obj.modifiedMask(obj.nsamples + (1:S)) = false;

            % 2. Spectral Filtering (Yule-Walker)
            % Filter data(:,P+1:end) — the lookahead-shifted window — for covariance
            % statistics only. This matches asr_process_version1: filter(B,A,data(:,range+P))
            % The carry prefix occupies indices 1:P, so P+1:end is the actual signal.
            obj.tProbe.start('spectral_shaping');
            [X_filt, obj.state.iir] = filter(obj.state.B, obj.state.A, data(:, P+1:end), obj.state.iir, 2);
            obj.tProbe.stop('spectral_shaping');

            % 3. Running Covariance Estimation
            obj.tProbe.start('moving_covariance');
            [Xcov, obj.state.cov] = moving_average(N, reshape(bsxfun(@times, ...
                reshape(X_filt, 1, C, []), reshape(X_filt, C, 1, [])), C*C, []), obj.state.cov);
            obj.tProbe.stop('moving_covariance');

            update_at = min(obj.stepsize:obj.stepsize:(size(Xcov,2)+obj.stepsize-1), size(Xcov,2));
            if isempty(obj.state.last_R)
                update_at          = [1, update_at];
                obj.state.last_R   = eye(C);
            end

            Xcov_mat = reshape(Xcov(:, update_at), C, C, []);
            last_n   = 0;

            % 4. Main Reconstruction & Adaptation Loop
            for j = 1:length(update_at)
                n = update_at(j);

                % Local covariance from unfiltered data window — same as vanillaASR.
                % Xcov_mat (filtered moving average) is only used for EMA adaptation,
                % NOT for the keep/trivial threshold comparison. Using the filtered
                % moving-average covariance for thresholding causes scale mismatch
                % with T (which is calibrated in the unfiltered signal space) and
                % results in every window being flagged as non-trivial.
                obj.tProbe.start('subspace_eig');
                % eig of the YW-filtered moving-average covariance — same spectral
                % space as T (which was calibrated on YW-filtered signal).
                % D and T*V are comparable because both live in the filtered space.
                [V, D_mat] = eig(Xcov_mat(:,:,j));
                [D, order] = sort(reshape(diag(D_mat), 1, C));
                V          = V(:, order);   % ascending order, col 1 = smallest

                % Local unfiltered cov for Riemannian telemetry only (different space is fine)
                idx_win = max(1, n-N+1)+P : n+P;
                Xw      = data(:, idx_win);
                Ct      = cov(Xw');
                Ct      = (Ct + Ct') / 2 + 1e-10 * eye(C);
                obj.tProbe.stop('subspace_eig');

                % Threshold comparison — both D and T*V in YW-filtered space ✓
                threshold_dir = sum((obj.T * V).^2);
                keep    = D < threshold_dir | (1:C) < (C - obj.maxdims * C);
                trivial = all(keep);

                % Snapshot M/T before any potential adaptation this step
                M_prev = obj.M;
                T_prev = obj.T;

                % FIX #6: removed j < length(update_at) so final window can adapt
                if trivial && j > 1 && mod(j, 20) == 0
                    subrange_idx = (last_n+1)+P : n+P;  % offset into carry-prefixed data

                    % FIX #5: approximateEntropy expects a vector — use mean channel
                    % power as a scalar proxy rather than passing the full matrix.
                    seg_power = mean(var(data(:, subrange_idx), 0, 2));
                    use_entropy_guard = (seg_power > 0);  % basic finite-signal check
                    if use_entropy_guard
                        % Only call approximateEntropy if it's on the path;
                        % fall back to always-adapt if the function is unavailable.
                        try
                            sig_scalar = mean(data(:, subrange_idx), 1); % [1 x T] mean
                            ent_ok = approximateEntropy(sig_scalar, 2) > 0.4;
                        catch
                            ent_ok = true;
                        end

                        if ent_ok
                            obj.tProbe.start('ema_adaptation');
                            obj.adaptBaseline(data(:, subrange_idx));
                            obj.tProbe.stop('ema_adaptation');

                            % Re-evaluate keep/trivial — V still from Xcov_mat (filtered space)
                            threshold_dir = sum((obj.T * V).^2);
                            keep    = D < threshold_dir | (1:C) < (C - obj.maxdims * C);
                            trivial = all(keep);
                        end
                    end
                end

                % Calculate Reconstruction Matrix R
                obj.tProbe.start('reconstruction_matrix');
                if ~trivial
                    R = real(obj.M * pinv(bsxfun(@times, keep', V' * obj.M)) * V');
                else
                    R = eye(C);
                end
                obj.tProbe.stop('reconstruction_matrix', struct('trivial', trivial));

                % Apply Blending (Raised-Cosine)
                if ~trivial || ~obj.state.last_trivial
                    % subrange in X_filt space (1..S); offset by P for data array
                    subrange      = (last_n+1):n;
                    subrange_data = subrange + P;
                    blend = (1 - cos(pi * (1:length(subrange)) / length(subrange))) / 2;
                    data(:, subrange_data) = ...
                        bsxfun(@times, blend,   R               * data(:, subrange_data)) + ...
                        bsxfun(@times, 1-blend, obj.state.last_R * data(:, subrange_data));

                    % Map back to original signal timeline
                    mask_idx = obj.nsamples + subrange;
                    mask_idx = mask_idx(mask_idx >= 1 & mask_idx <= obj.nsamples + S);
                    obj.modifiedMask(mask_idx) = ~trivial;
                end

                % Telemetry — Riemannian distance always measured from frozen
                % calibration baseline (obj.state.M), not the EMA-adapted obj.M.
                % FIX #7: this is intentional — kept as obj.state.M per design decision.
                delta_anchor = (obj.state.M \ Ct) / obj.state.M';
                delta_anchor = (delta_anchor + delta_anchor') / 2;
                riem_dist    = sqrt(sum(log(max(eig(delta_anchor), 1e-12)).^2));

                bufFill = min(1, size(obj.cleanBuffer,2) / ...
                    max(1, round(obj.cleanBuffLength * obj.srate)));
                obj.emaProbe.update(Ct, V, D, R, trivial, ...
                                    riem_dist, obj.M, obj.T, M_prev, T_prev, bufFill);

                % State Update
                last_n                 = n;
                obj.state.last_R       = R;
                obj.state.last_trivial = trivial;
            end

            % 5. Output Management
            % Blending was applied to data(:,subrange_data) = data(:,subrange+P)
            % so carry region data(:,1:P) is untouched boundary reflection.
            % Strip carry prefix only → exactly S samples out per chunk.
            % Lookahead is implicit: we output the current chunk delayed by the
            % boundary carry (same convention as asr_process_version1).
            obj.state.carry = data(:, (end-P+1):end);
            cleanedSignal   = data(:, P+1:end);
            obj.nsamples    = obj.nsamples + S;
        end

    end % end public methods

    methods (Access = private)

        function adaptBaseline(obj, segment)
            % adaptBaseline - Robust L1-Median Evolution with GPU/CPU dispatch.

            % 1. Rolling Buffer Management
            max_samples     = round(obj.cleanBuffLength * obj.srate);
            obj.cleanBuffer = [obj.cleanBuffer, segment];
            if size(obj.cleanBuffer, 2) > max_samples
                obj.cleanBuffer = obj.cleanBuffer(:, end-max_samples+1:end);
            end

            % 2. Block Covariance Accumulation (GPU if available, else CPU)
            % FIX #8: canUseGPU() guard prevents hard-error on CPU-only machines
            % FIX #3/#4: dimension mismatch and normalisation corrected —
            %   accumulate all outer products, divide by actual sample count.
            Xbuf     = obj.cleanBuffer';          % [N_samples x nchans]
            N_samples = size(Xbuf, 1);

            use_gpu = canUseGPU();
            if use_gpu
                Xbuf = gpuArray(Xbuf);
            end

            % Accumulate [nchans x nchans] covariance blocks and track sample count
            U_acc  = zeros(obj.nchans, obj.nchans, 'like', Xbuf);
            n_used = 0;

            for k = 1:obj.blocksize
                idx = k:obj.blocksize:N_samples;
                if isempty(idx), continue; end
                Xk     = Xbuf(idx, :);               % [n_k x nchans]
                U_acc  = U_acc + (Xk' * Xk);         % [nchans x nchans]
                n_used = n_used + size(Xk, 1);
            end

            if n_used == 0, return; end

            % Reshape each block's normalised outer product for geometric median
            % Build [n_blocks x nchans^2] matrix expected by block_geometric_median
            n_blocks  = ceil(N_samples / obj.blocksize);
            U_blocks  = zeros(n_blocks, obj.nchans^2, 'like', Xbuf);
            b = 0;
            for k = 1:obj.blocksize
                idx = k:obj.blocksize:N_samples;
                if isempty(idx), continue; end
                b = b + 1;
                Xk           = Xbuf(idx, :);
                cov_k        = (Xk' * Xk) / size(Xk, 1);  % FIX #4: divide by actual count
                U_blocks(b,:) = cov_k(:)';
            end
            U_blocks = U_blocks(1:b, :);   % trim unused rows

            % 3. Geometric Median → new Mixing Matrix
            if use_gpu
                U_med = gather(block_geometric_median(U_blocks));
            else
                U_med = block_geometric_median(U_blocks);
            end
            M_new = sqrtm(real(reshape(U_med, obj.nchans, obj.nchans)));

            % 4. EMA Updates
            obj.M = (1 - obj.beta1) * obj.M + obj.beta1 * M_new;

            % Recalculate T based on the adapted M and the current clean buffer
            T_new = emaASR.computeThreshold(obj.window_len, obj.srate, ...
                size(obj.cleanBuffer, 2), obj.M, obj.cleanBuffer', obj.nchans, ...
                [], [], [], obj.cutoff);

            obj.T = (1 - obj.beta2) * obj.T + obj.beta2 * T_new;
        end

    end % end private methods

    % ================================================================
    % Private static helpers — inlined calibration pipeline
    % (replaces: asr_calibrate_Version1, applytheoriginalYuleWalker,
    %  applyYuleWalkFilter, getMixingMatrixForReconstruction,
    %  calculatePerComponentThreshold, intializeCalibrationStates)
    % ================================================================
    methods (Static, Access = private)

        function state = calibrateBaseline(X, srate, cutoff, blocksize, window_len)
            % CALIBRATEBASELINE - Full calibration pipeline, self-contained.
            %   Inlines asr_calibrate_Version1 with no external file dependencies.
            %
            %   Steps:
            %     1. Auto-scale blocksize to available memory
            %     2. Spectral shaping (Yule-Walker IIR)
            %     3. Block geometric median → mixing matrix M
            %     4. Per-component threshold matrix T
            %     5. Build and return calibration state struct

            [C, S] = size(X);

            % Auto-scale blocksize to available memory
            try
                maxmem    = hlp_memfree / (2^21);
                blocksize = max(blocksize, ceil((C*C*S*8*3*2) / (maxmem*(2^21))));
            catch
                % hlp_memfree not available — use blocksize as given
            end

            % Step 1: Spectral shaping (Yule-Walker IIR)
            [B, A]           = emaASR.buildYuleWalkerFilter(srate);
            X(~isfinite(X))  = 0;
            [X, iirstate]    = filter(B, A, double(X), [], 2);

            % Step 2: Block geometric median covariance → mixing matrix M
            [X, M] = emaASR.computeMixingMatrix(X, blocksize, S, C);

            % Step 3: Per-component threshold matrix T
            T = emaASR.computeThreshold(window_len, srate, S, M, X, C, ...
                [], [], [], cutoff);

            % Step 4: Build state struct
            state = emaASR.buildCalibrationState(M, T, B, A, iirstate, X);
        end

        % ----------------------------------------------------------------
        function [B, A] = buildYuleWalkerFilter(srate)
            % BUILDYULEWALKERFILTER - Spectral shaping IIR filter.
            %   Inlines applyYuleWalkFilter + applytheoriginalYuleWalker.
            %   Tries yulewalk (Signal Processing Toolbox); falls back to
            %   pre-computed coefficients for standard sampling rates.
            try
                freqvals = [[0 2 3 13 16 40 min(80, srate/2-1)] * 2/srate, 1];
                amps     = [3 0.75 0.33 0.33 1 1 3 3];
                if srate < 80
                    freqvals(end-2) = [];
                    amps(end)       = [];
                end
                [B, A] = yulewalk(8, freqvals, amps);
            catch
                % Pre-computed coefficients for standard sampling rates
                switch srate
                    case 100
                        B = [0.9314233528641650 -1.0023683814963549 -0.4125359862018213  0.7631567476327510  0.4160430392910331 -0.6549131038692215 -0.0372583518046807  0.1916268458752655  0.0462411971592346];
                        A = [1.0000000000000000 -0.4544220180303844 -1.0007038682936749  0.5374925521337940  0.4905013360991340 -0.4861062879351137 -0.1995986490699414  0.1830048420730026  0.0457678549234644];
                    case 128
                        B = [1.1027301639165037 -2.0025621813611867  0.8942119516481342  0.1549979524226999  0.0192366904488084  0.1782897770278735 -0.5280306696498717  0.2913540603407520 -0.0262209802526358];
                        A = [1.0000000000000000 -1.1042042046423233 -0.3319558528606542  0.5802946221107337 -0.0010360013915635  0.0382167091925086 -0.2609928034425362  0.0298719057761086  0.0935044692959187];
                    case 200
                        B = [1.4489483325802353 -2.6692514764802775  2.0813970620731115 -0.9736678877049534  0.1054605060352928 -0.1889101692314626  0.6111331636592364 -0.3616483013075088  0.1834313060776763];
                        A = [1.0000000000000000 -0.9913236099393967  0.3159563145469344 -0.0708347481677557 -0.0558793822071149 -0.2539619026478943  0.2473056615251193 -0.0420478437473110  0.0077455718334464];
                    case 256
                        B = [1.7587013141770287 -4.3267624394458641  5.7999880031015953 -6.2396625463547508  5.3768079046882207 -3.7938218893374835  2.1649108095226470 -0.8591392569863763  0.2569361125627988];
                        A = [1.0000000000000000 -1.7008039639301735  1.9232830391058724 -2.0826929726929797  1.5982638742557307 -1.0735854183930011  0.5679719225652651 -0.1886181499768189  0.0572954115997261];
                    case 300
                        B = [1.9153920676433143  -5.7748421104926795   9.1864764859103936 -10.7350356619363630   9.6423672437729007  -6.6181939699544277   3.4219421494177711  -1.2622976569994351   0.2968423019363821];
                        A = [1.0000000000000000 -2.3143703322055491  3.2222567327379434 -3.6030527704320621  2.9645154844073698 -1.8842615840684735  0.9222455868758080 -0.3103251703648485  0.0634586449896364];
                    case 500
                        B = [2.3133520086975823 -11.9471223009159130  29.1067166493384340 -43.7550171007238190  44.3385767452216370 -30.9965523846388000  14.6209883020737190  -4.2743412400311449   0.5982553583777899];
                        A = [1.0000000000000000  -4.6893329084452580  10.5989986701080210 -14.9691518101365230  14.3320358399731820  -9.4924317069169977   4.2425899618982656  -1.1715600975178280   0.1538048427717476];
                    case 512
                        B = [2.3275475636130865 -12.2166478485960430  30.1632789058248850 -45.8009842020820410  46.7261263011068880 -32.7796858196767220  15.4623349612560630  -4.5019779685307473   0.6242733481676324];
                        A = [1.0000000000000000  -4.7827378944258703  10.9780696236622980 -15.6795187888195360  15.1281978667576310 -10.0632079834518220   4.5014690636505614  -1.2394100873286753   0.1614727510688058];
                    otherwise
                        error('emaASR:noFilter', ...
                            ['No pre-computed Yule-Walker filter for srate=%d Hz. ' ...
                             'Install the Signal Processing Toolbox (yulewalk) or ' ...
                             'resample to: 100, 128, 200, 256, 300, 500, 512 Hz.'], srate);
                end
            end
        end

        % ----------------------------------------------------------------
        function [Xt, M] = computeMixingMatrix(X, blocksize, S, C)
            % COMPUTEMIXINGMATRIX - Block geometric median covariance → sqrtm.
            %   Inlines getMixingMatrixForReconstruction.
            %   X input  : [C x S] (filter output, rows = channels)
            %   Xt output: [S x C] (transposed, rows = samples)
            Xt = X';   % [S x C]
            U  = zeros(length(1:blocksize:S), C*C);
            for k = 1:blocksize
                range = min(S, k:blocksize:(S+k-1));
                U = U + reshape(bsxfun(@times, ...
                    reshape(Xt(range,:), [], 1, C), ...
                    reshape(Xt(range,:), [], C, 1)), size(U));
            end
            M = sqrtm(real(reshape(block_geometric_median(U / blocksize), C, C)));
        end

        % ----------------------------------------------------------------
        function T = computeThreshold(window_len, srate, S, M, X, C, ...
                window_overlap, min_clean_fraction, max_dropout_fraction, cutoff)
            % COMPUTETHRESHOLD - Per-component ASR threshold matrix.
            %   Inlines calculatePerComponentThreshold.
            %   X : [S x C] (rows = samples)
            if nargin < 7  || isempty(window_overlap),       window_overlap       = 0.66; end
            if nargin < 8  || isempty(min_clean_fraction),   min_clean_fraction   = 0.25; end
            if nargin < 9  || isempty(max_dropout_fraction), max_dropout_fraction = 0.1;  end
            if nargin < 10 || isempty(cutoff),               cutoff               = 20;   end

            N = round(window_len * srate);
            if S < N
                error('emaASR:tooShort', ...
                    'Calibration data too short for window_len=%.2f s at %d Hz.', ...
                    window_len, srate);
            end

            [V, ~] = eig(M);
            Xv     = abs(X * V);   % [S x C]

            mu  = zeros(1, C);
            sig = zeros(1, C);
            for c = C:-1:1
                rms_sq = Xv(:, c) .^ 2;
                wins   = bsxfun(@plus, ...
                    round(1 : N*(1-window_overlap) : S-N)', (0:N-1));
                rmsw   = sqrt(sum(rms_sq(wins), 2) / N);
                [mu(c), sig(c)] = fit_eeg_distribution( ...
                    rmsw, min_clean_fraction, max_dropout_fraction);
            end
            T = diag(mu + cutoff * sig) * V';
        end

        % ----------------------------------------------------------------
        function state = buildCalibrationState(M, T, B, A, iirstate, X)
            % BUILDCALIBRATIONSTATE - Pack calibration outputs into state struct.
            %   Inlines intializeCalibrationStates.
            %   X : [S x C] (rows = samples, from computeMixingMatrix output)
            state = struct( ...
                'M',            M,           ...
                'T',            T,           ...
                'B',            B,           ...
                'A',            A,           ...
                'iir',          iirstate,    ...
                'cov',          [],          ...
                'carry',        [],          ...
                'last_R',       [],          ...
                'last_trivial', true,        ...
                'last_Xcov',    [],          ...
                'cleanBuffer',  X);

            % Entropy gating bookkeeping
            state.entropy_thr  = 0.4;
            state.entropy_skip = 20;
            state.entropy_cnt  = 0;
            state.entropy_low  = 0.4;
            state.entropy_high = 1.2;

            % PCA bookkeeping
            state.pca_skip = 5;
            state.pca_cnt  = 0;

            % Eigenspace bookkeeping
            state.eigen_threshold = 0.05;
            state.last_V          = [];
            state.last_lam        = [];
        end

    end % static private methods

end
