classdef hmoASR < handle
    % hmoASR - Hierarchical Model-Order Artifact Subspace Reconstruction.
    %   Online EEG artifact removal with a continuously-updated clean covariance
    %   model. Each processing window is gated as clean or artefactual by comparing
    %   its per-PC RMS against a calibrated threshold. Clean windows update the
    %   global covariance model (Welford online) and eigenbasis (early-eigen reuse);
    %   artefactual windows are reconstructed via an oblique projector in the clean
    %   subspace.
    %
    %   Calibration gate : whitened covariance deviation (z_dev) — Frobenius norm
    %                      of (sqrtC \ Cw / sqrtC - I) z-scored online.
    %                      Calibration data is NOT Level-1 pruned (Wpp = W).
    %   Processing gate  : Level-1 PCA z-prune + channel-RMS z-score.
    %   Output convention: P+1:end  (consistent with emaASR — no lookahead delay).
    %
    %   Reference:
    %       calibrate_hmo_asr.m / process_hmo_asr.m (Hübner et al., 2023)
    %
    %   Constructor signature (matches Experimenting.m harness):
    %       obj = hmoASR(srate, params)
    %       obj = hmoASR(srate, params, tProbe)

    % ================================================================
    properties (Access = public)

        % --- Core ---
        nchans
        srate
        nsamples    = 0
        modifiedMask

        % --- Algorithm Parameters ---
        cutoff              = 10        % SD multiplier for per-PC amplitude threshold
        window_len          = 0.5       % Processing window length (s)
        window_overlap      = 0.0       % Calibration window overlap fraction
        eigen_threshold     = 0.002     % Off-diagonal ratio triggering eigenbasis refresh
        z_clean_range       = [-3 5]    % z-score gate [lo hi]
        maxdims             = 1         % Max PCs to reject per window (int or fraction <1)
        lookahead           = 0.25      % Lookahead buffer (s)
        stepsize            = 0.5       % Hop between processing windows (s)

        % --- Spectral Shaping Filter (Yule-Walker) ---
        B                               % Numerator coefficients
        A                               % Denominator coefficients

        % --- Clean Model State ---
        n_tot   = 0                     % Effective clean sample count
        m_ch                            % Running channel mean  [C x 1]
        C_ch                            % Running clean covariance  [C x C]
        E                               % Current eigenbasis  [C x C]
        sqrtC                           % sqrtm(C_ch) in current basis  [C x C]
        mu_pc                           % Per-PC RMS mean  [C x 1]
        var_pc                          % Per-PC RMS variance  [C x 1]
        Tpc                             % Per-PC amplitude thresholds  [C x 1]

        % --- Channel-RMS Gate State (initialised in calibrate, updated in process) ---
        n_rms   = 0
        ch_mu
        ch_var

        % --- Runtime Buffers ---
        carry   = []                    % Lookahead boundary buffer  [C x P]
        iir_state = []                  % IIR spectral filter state

        % --- Telemetry ---
        tProbe                          % timeProbe instance
        probe                           % hmoASRprobe instance

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = hmoASR(srate, params, tProbe)
            % HMOASR
            %   srate  : Sampling rate in Hz (required)
            %   params : Struct of parameter overrides (optional)
            %   tProbe : timeProbe instance (optional)

            if nargin < 1 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('hmoASR:invalidSrate', 'srate must be a positive numeric scalar.');
            end
            obj.srate = srate;

            if nargin >= 3 && ~isempty(tProbe)
                obj.tProbe = tProbe;
            else
                obj.tProbe = timeProbe();
            end

            if nargin >= 2 && isstruct(params)
                f = fieldnames(params);
                for i = 1:numel(f)
                    if isprop(obj, f{i})
                        obj.(f{i}) = params.(f{i});
                    end
                end
            end

            [obj.B, obj.A] = hmoASR.buildSpectralFilter(srate);
            obj.modifiedMask = false(1, 0);
        end

        % ── Calibration ──────────────────────────────────────────────
        function calibrate(obj, X)
            % CALIBRATE - Establish initial clean model from resting EEG.
            %
            %   Faithful to calibrate_hmo_asr.m:
            %     - No Level-1 prune during calibration (Wpp = W)
            %     - Whitened covariance deviation (z_dev) acceptance gate:
            %         Cw    = cov(W, 1)
            %         Delta = sqrtC \ Cw / sqrtC
            %         dev   = ||Delta - I||_F
            %         Accept window if z_dev in z_clean_range
            %     - Welford parallel mean/cov update on accepted windows
            %     - Early-eigen reuse for basis refresh
            %     - Per-PC RMS stats → Tpc threshold

            if isempty(X) || ~isnumeric(X)
                error('hmoASR:calibrate', 'Calibration data must be non-empty numeric.');
            end
            X = obj.ensureOrientation(X);
            fprintf('--- [hmoASR] Calibrating ---\n');
            obj.tProbe.start('calibration');

            [C, S] = size(X);
            obj.nchans = C;
            X(~isfinite(X)) = 0;

            % Spectral shaping — preserve IIR state for seamless process() startup
            [Xf, iirstate_cal] = filter(obj.B, obj.A, double(X), [], 2);
            Xf = Xf.';   % [S x C]

            Nw  = max(1, round(obj.window_len * obj.srate));
            hop = max(1, round(Nw * (1 - obj.window_overlap)));
            win_starts = 1 : hop : max(1, S - Nw + 1);

            % Running clean model accumulators
            n_tot = 0;  m_ch = zeros(C,1);  C_ch = zeros(C);
            E_prev = [];  sqrtC_run = eye(C);  mu_pc = [];  var_pc = [];

            % Whitened-deviation gate accumulators
            n_dev = 0;  mu_dev = 0;  var_dev = 1;

            for ws = 1:numel(win_starts)
                idx = win_starts(ws) : win_starts(ws) + Nw - 1;
                if idx(end) > S, break; end

                W = Xf(idx, :);   % [Nw x C] — no Level-1 prune during calibration

                % --- Ensure sqrtC exists before z_dev gate can fire ---
                if isempty(E_prev)
                    % Bootstrap: first window always accepted (no baseline yet)
                    n_star = size(W, 1);
                    if n_star < 2, continue; end
                    m_star = mean(W, 1).';
                    Xc     = W - m_star.';
                    C_ch   = (Xc.' * Xc) / (n_star - 1);
                    m_ch   = m_star;  n_tot = n_star;

                    [E_prev, Dc] = eig((C_ch + C_ch.') / 2 + 1e-10*eye(C));
                    lam          = max(0, real(diag(Dc)));
                    sqrtC_run    = E_prev * diag(sqrt(lam)) * E_prev.';

                    Y      = (E_prev.' * W.').';
                    pc_rms = sqrt(mean(Y.^2, 1)).';
                    mu_pc  = pc_rms;  var_pc = ones(C,1) * 1e-6;
                    continue;
                end

                % --- Whitened covariance deviation gate ---
                Cw    = cov(W, 1);
                Delta = sqrtC_run \ Cw / sqrtC_run;
                dev   = norm(Delta - eye(C), 'fro');

                if n_dev == 0
                    mu_dev = dev;  var_dev = 1e-6;  n_dev = 1;
                    z_dev  = 0;
                else
                    z_dev = (dev - mu_dev) / sqrt(var_dev + 1e-12);
                end

                if z_dev < obj.z_clean_range(1) || z_dev > obj.z_clean_range(2)
                    continue;
                end
                [mu_dev, var_dev, n_dev] = hmoASR.updateMeanVar(mu_dev, var_dev, n_dev, dev);

                % --- Welford update ---
                n_star = size(W, 1);
                if n_star < 2, continue; end
                m_star = mean(W, 1).';
                Xc     = W - m_star.';
                C_star = (Xc.' * Xc) / (n_star - 1);

                n_new = n_tot + n_star;
                delta = m_ch - m_star;
                m_new = m_ch + (m_star - m_ch) * (n_star / n_new);
                C_new = (((n_tot-1)*C_ch) + ((n_star-1)*C_star) + ...
                         (n_tot*n_star/n_new) * (delta*delta.')) / (n_new-1);
                n_tot = n_new;  m_ch = m_new;  C_ch = C_new;

                % --- Early-eigen reuse ---
                [E_prev, lam, mu_pc, var_pc] = hmoASR.refreshEigenbasis( ...
                    E_prev, C_ch, obj.eigen_threshold, mu_pc, var_pc);
                sqrtC_run = E_prev * diag(sqrt(lam)) * E_prev.';

                % --- Per-PC RMS stats ---
                Y      = (E_prev.' * W.').';
                pc_rms = sqrt(mean(Y.^2, 1)).';
                if isempty(mu_pc)
                    mu_pc = pc_rms;  var_pc = ones(C,1) * 1e-6;
                end
                [mu_pc, var_pc, ~] = hmoASR.updateMeanVar(mu_pc, var_pc, ws, pc_rms, true);
            end

            % Fallbacks if no clean windows were accepted
            if isempty(E_prev) || n_tot == 0
                E_prev    = eye(C);
                lam       = max(1e-12, real(diag((Xf.' * Xf) / max(1, S))));
                sqrtC_run = E_prev * diag(sqrt(lam)) * E_prev.';
            end
            if isempty(mu_pc)
                lam_now = max(0, real(diag(E_prev.' * C_ch * E_prev)));
                mu_pc   = sqrt(lam_now);  var_pc = ones(C,1) * 1e-6;
            end

            % Store clean model
            obj.n_tot  = n_tot;   obj.m_ch  = m_ch;   obj.C_ch  = C_ch;
            obj.E      = E_prev;  obj.sqrtC = sqrtC_run;
            obj.mu_pc  = mu_pc;   obj.var_pc = var_pc;
            obj.Tpc    = mu_pc + obj.cutoff * sqrt(max(var_pc, 0));

            % Channel-RMS gate initial state for process()
            obj.n_rms  = 0;  obj.ch_mu = zeros(C,1);  obj.ch_var = ones(C,1);

            % Preserve IIR state — avoids startup transient on first process() call
            obj.iir_state = iirstate_cal;
            obj.carry     = [];

            obj.probe = hmoASRprobe(C);
            obj.probe.setBaseline(obj);

            obj.nsamples     = 0;
            obj.modifiedMask = false(1, 0);

            obj.tProbe.stop('calibration');
            fprintf('--- [hmoASR] Calibration complete (n_tot=%d) ---\n', n_tot);
        end

        % ── Process ──────────────────────────────────────────────────
        function cleanedSignal = process(obj, data)
            % PROCESS - Online artifact rejection with adaptive clean model.
            %   Faithful to process_hmo_asr.m.

            if isempty(obj.E)
                error('hmoASR:notCalibrated', 'Run calibrate() before process().');
            end
            if isempty(data) || ~isnumeric(data)
                error('hmoASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            data = obj.ensureOrientation(data);

            obj.tProbe.start('process');
            [C, S] = size(data);

            % Convert fraction maxdims → integer count
            mdims = obj.maxdims;
            if mdims > 0 && mdims < 1, mdims = round(C * mdims); end

            Nw  = max(1, round(obj.window_len * obj.srate));
            hop = max(1, round(obj.stepsize   * obj.srate));
            P   = max(0, round(obj.lookahead  * obj.srate));

            obj.modifiedMask(obj.nsamples + (1:S)) = false;

            % Lookahead boundary (linear extrapolation seed — matches process_hmo_asr.m)
            if isempty(obj.carry)
                if S >= 2
                    obj.carry = repmat(2*data(:,1) - data(:,2), 1, P);
                else
                    obj.carry = repmat(data(:,1), 1, P);
                end
            end
            data_ext = [obj.carry, double(data)];
            data_ext(~isfinite(data_ext)) = 0;

            obj.tProbe.start('spectral_shaping');
            [Xf, obj.iir_state] = filter(obj.B, obj.A, data_ext, obj.iir_state, 2);
            Xf = Xf.';   % [P+S x C]
            obj.tProbe.stop('spectral_shaping');

            totalS     = size(Xf, 1);
            win_starts = (P+1) : hop : (totalS - Nw + 1);
            cleaned    = data_ext;

            for w = 1:numel(win_starts)
                a = win_starts(w);
                b = a + Nw - 1;

                W_raw  = data_ext(:, a:b).';
                W_filt = Xf(a:b, :);

                % Level-1 PCA z-prune
                obj.tProbe.start('gating');
                [Wpp, ~] = hmoASR.pcaEigZPrune(W_filt, 1.5);

                % Level-2 channel-RMS z-gate
                ch_rms = sqrt(mean(Wpp.^2, 1)).';
                [obj.ch_mu, obj.ch_var, obj.n_rms] = ...
                    hmoASR.updateMeanVar(obj.ch_mu, obj.ch_var, obj.n_rms, ch_rms);
                z_ch     = (ch_rms - obj.ch_mu) ./ max(1e-12, sqrt(obj.ch_var));
                is_clean = all(z_ch >= obj.z_clean_range(1) & z_ch <= obj.z_clean_range(2));
                obj.tProbe.stop('gating');

                % Clean-model update (Eq. 4-6)
                if is_clean
                    obj.tProbe.start('model_update');
                    obj.updateCleanModel(Wpp, w);
                    obj.tProbe.stop('model_update');
                end

                % Reconstruction (Eq. 10-11)
                obj.tProbe.start('reconstruction_matrix');
                [Rproj, Ek, lam, Tproj, reject] = obj.buildProjector(W_filt, Nw, mdims);
                obj.tProbe.stop('reconstruction_matrix');

                Wclean  = (Rproj * W_raw.').';
                trivial = ~any(reject);

                cleaned(:, a:b) = Wclean.';
                if ~trivial
                    mask_idx = obj.nsamples + (a:b) - P;
                    mask_idx = mask_idx(mask_idx >= 1 & mask_idx <= obj.nsamples + S);
                    obj.modifiedMask(mask_idx) = true;
                end

                riem_dist = obj.calcGeodesic(obj.C_ch, (W_filt.'*W_filt)/max(1,Nw));
                obj.probe.update(obj, (W_filt.'*W_filt)/max(1,Nw), ...
                    Ek, lam, Tproj, Rproj, reject, riem_dist);
            end

            % Carry = last P samples of this chunk (prepended to next call)
            obj.carry     = data_ext(:, end-P+1:end);
            cleanedSignal = cleaned(:, P+1:end);   % strip carry → S samples, no delay
            obj.nsamples  = obj.nsamples + S;

            obj.tProbe.stop('process');
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function X = ensureOrientation(obj, X)
            if ~isnumeric(X) || isempty(X)
                error('hmoASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            if ~isempty(obj.nchans)
                if size(X,1) ~= obj.nchans && size(X,2) == obj.nchans
                    X = X.';
                elseif size(X,1) ~= obj.nchans
                    error('hmoASR:channelMismatch', ...
                        'Expected %d channels, got [%d x %d].', ...
                        obj.nchans, size(X,1), size(X,2));
                end
            else
                if size(X,1) > size(X,2), X = X.'; end
            end
        end

        function updateCleanModel(obj, Wpp, win_idx)
            % UPDATECLEANMODEL - Welford update + early-eigen reuse (Eq. 4-6).
            n_star = size(Wpp, 1);
            if n_star < 2, return; end

            m_star = mean(Wpp, 1).';
            Xc     = Wpp - m_star.';
            C_star = (Xc.' * Xc) / (n_star - 1);

            n_new = obj.n_tot + n_star;
            if obj.n_tot == 0
                m_new = m_star;  C_new = C_star;
            else
                delta = obj.m_ch - m_star;
                m_new = obj.m_ch + (m_star - obj.m_ch) * (n_star / n_new);
                C_new = (((obj.n_tot-1)*obj.C_ch) + ((n_star-1)*C_star) + ...
                         (obj.n_tot*n_star/n_new) * (delta*delta.')) / (n_new-1);
            end
            obj.n_tot = n_new;  obj.m_ch = m_new;  obj.C_ch = C_new;

            [obj.E, lam, obj.mu_pc, obj.var_pc] = hmoASR.refreshEigenbasis( ...
                obj.E, obj.C_ch, obj.eigen_threshold, obj.mu_pc, obj.var_pc);
            obj.sqrtC = obj.E * diag(sqrt(lam)) * obj.E.';

            Y      = (obj.E.' * Wpp.').';
            pc_rms = sqrt(mean(Y.^2, 1)).';
            if isempty(obj.mu_pc)
                obj.mu_pc  = pc_rms;
                obj.var_pc = ones(obj.nchans,1) * 1e-6;
            end
            [obj.mu_pc, obj.var_pc, ~] = hmoASR.updateMeanVar( ...
                obj.mu_pc, obj.var_pc, win_idx, pc_rms, true);
            obj.Tpc = obj.mu_pc + obj.cutoff * sqrt(max(obj.var_pc, 0));
        end

        function [Rproj, Ek, lam, Tproj, reject] = buildProjector(obj, W_filt, Nw, mdims)
            % BUILDPROJECTOR - Oblique projector in clean subspace (Eq. 10-11).
            COVw     = (W_filt.' * W_filt) / max(1, Nw);
            [Ek, Dk] = eig((COVw + COVw.') / 2);
            lam      = real(diag(Dk));

            if isempty(obj.E)
                Tproj = obj.Tpc .^ 2;
            else
                Q     = obj.E.' * Ek;
                Tproj = (Q.^2).' * (obj.Tpc .^ 2);
            end
            reject = lam > Tproj;

            C = obj.nchans;
            if mdims < C
                [~, ord]   = sort(lam, 'ascend');
                keep_force = ord(1:(C - mdims));
                reject(keep_force) = false;
            end

            G     = Ek.' * obj.sqrtC;
            Gt    = G;  Gt(reject,:) = 0;
            Rproj = obj.sqrtC * pinv(Gt) * Ek.';
        end

        function dist = calcGeodesic(~, C1, C2)
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

    % ================================================================
    methods (Static, Access = private)

        function [B, A] = buildSpectralFilter(srate)
            freqvals = [[0 2 3 13 16 40 min(80, srate/2-1)] * 2/srate, 1];
            amps     = [3 0.75 0.33 0.33 1 1 3 3];
            if srate < 80, freqvals(end-2) = [];  amps(end) = []; end
            try
                [B, A] = yulewalk(8, freqvals, amps);
            catch
                [B, A] = butter(4, 2*[0.5 min(45, srate/2-2)] / srate, 'bandpass');
            end
        end

        function [E_new, lam, mu_pc, var_pc] = refreshEigenbasis(E_prev, C_ch, thr, mu_pc, var_pc)
            % REFRESHEIGENBASIS - Early-eigen reuse: only refactor if basis has drifted.
            if isempty(E_prev)
                need = true;
            else
                Rtest = E_prev.' * C_ch * E_prev;
                offd  = sum(abs(Rtest(:))) - sum(abs(diag(Rtest)));
                diagv = sum(abs(diag(Rtest))) + 1e-12;
                need  = (offd / diagv) >= thr;
            end
            if need
                E_old = E_prev;
                if isempty(E_old), E_old = eye(size(C_ch,1)); end
                [E_new, Dc] = eig((C_ch + C_ch.')/2);
                lam = max(0, real(diag(Dc)));
                if ~isempty(mu_pc)
                    [mu_pc, var_pc] = hmoASR.projectStatsNewBasis(mu_pc, var_pc, E_old, E_new);
                end
            else
                E_new = E_prev;
                lam   = max(0, real(diag(E_new.' * C_ch * E_new)));
            end
        end

        function [Wclean, kept] = pcaEigZPrune(W, zcut)
            COV  = (W.' * W) / max(1, size(W,1));
            [E, D] = eig((COV + COV.')/2);
            lam  = real(diag(D));
            z    = (lam - mean(lam)) / (std(lam) + 1e-12);
            kept = z <= zcut;
            E2   = E(:, kept);
            if isempty(E2)
                [~,i] = min(lam);  E2 = E(:,i);
                kept = false(size(kept));  kept(i) = true;
            end
            Wclean = (E2 * (E2.' * W.')).';
        end

        function [mu, varv, n] = updateMeanVar(mu, varv, n, x, robust_init)
            if nargin < 5, robust_init = false; end
            x = x(:);
            if n == 0
                mu = x;  varv = ones(size(x)) * (robust_init * 1e-6);  n = 1;  return;
            end
            n_new  = n + 1;
            delta  = x - mu;
            mu_new = mu + delta ./ n_new;
            M2_old = varv .* max(1, n-1);
            M2_new = M2_old + delta .* (x - mu_new);
            varv   = max(M2_new ./ max(1, n_new-1), 0);
            mu = mu_new;  n = n_new;
        end

        function [mu_new, var_new] = projectStatsNewBasis(mu_old, var_old, E_old, E_new)
            Q       = E_old.' * E_new;
            s2_old  = mu_old .^ 2 + var_old;
            s2_new  = (Q.^2).' * s2_old;
            var_new = (Q.^2).' * var_old;
            mu_new  = sqrt(max(s2_new - var_new, 0));
        end

    end % static private methods

end
