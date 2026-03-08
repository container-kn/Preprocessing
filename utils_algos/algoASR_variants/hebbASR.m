classdef hebbASR < handle
    % hebbASR - Hebbian/Anti-Hebbian Artifact Subspace Reconstruction (PSP/PSW).
    %   Extends hmoASR calibration with online Hebbian weight learning. The
    %   eigenbasis is no longer recomputed via explicit eig() each window; instead
    %   feedforward weights W and lateral weights M are updated via the PSP or PSW
    %   gradient rules (Tsai et al., 2024), and E is derived from QR(M\W)'.
    %
    %   Two learning rules:
    %     'psp' — Pehlevan-Sreekumar-Pfister (Algorithm 1):
    %               Y   = M \ (W * Xk)
    %               W  ← (1-rW)·W  + rW·(Y·Xk')
    %               M  ← (1-rM)·M  + rM·(Y·Y')
    %               rW = 2η,  rM = rW/(2τ)
    %     'psw' — PSW with whitening constraint (Algorithm 2):
    %               Y   = M \ (W * Xk)
    %               W  ← (1-rW)·W  + rW·(Y·Xk')
    %               M  ← (M - rM·I) + rM·(Y·Y')
    %               rW = η,  rM = rW/τ
    %
    %   Calibration: full hmoASR-style HMO calibration (with z_dev gate),
    %                then seeds W = E.' and M = I.
    %   Process:     same gating + Welford update as hmoASR, but eigenbasis
    %                derived from Hebbian weights rather than explicit eig.
    %   Output:      P+1:end (consistent with emaASR/hmoASR).
    %
    %   Reference:
    %       calibrate_hebb_asr.m / process_hebb_asr.m (Tsai et al., 2024)
    %
    %   Constructor signature (matches Experimenting.m harness):
    %       obj = hebbASR(srate, params)
    %       obj = hebbASR(srate, params, tProbe)

    % ================================================================
    properties (Access = public)

        % --- Core ---
        nchans
        srate
        nsamples    = 0
        modifiedMask

        % --- Algorithm Parameters (shared with hmoASR) ---
        cutoff              = 10
        window_len          = 0.5
        window_overlap      = 0.0
        eigen_threshold     = 0.002
        z_clean_range       = [-3 5]
        maxdims             = 1
        lookahead           = 0.25
        stepsize            = 0.5

        % --- Hebbian Hyperparameters ---
        hebb_mode   = 'psw'   % 'psp' or 'psw'
        eta         = 0.05    % Base learning rate (η)
        tau         = 10      % Lateral rate ratio (τ)
        rW                    % Derived feedforward learning rate
        rM                    % Derived lateral learning rate
        hebb_step   = 0       % Cumulative Hebbian update counter

        % --- Hebbian Weight Matrices ---
        W                     % Feedforward weights  [C x C]
        M                     % Lateral weights      [C x C]

        % --- Spectral Shaping Filter ---
        B
        A

        % --- Clean Model State (identical to hmoASR) ---
        n_tot   = 0
        m_ch
        C_ch
        E                     % Current eigenbasis (derived from W, M each step)
        sqrtC
        mu_pc
        var_pc
        Tpc

        % --- Channel-RMS Gate State ---
        n_rms   = 0
        ch_mu
        ch_var

        % --- Runtime Buffers ---
        carry     = []
        iir_state = []

        % --- Telemetry ---
        tProbe
        probe                 % hebbASRprobe instance

    end % properties

    % ================================================================
    methods

        % ── Constructor ──────────────────────────────────────────────
        function obj = hebbASR(srate, params, tProbe)
            % HEBBASR
            %   srate  : Sampling rate in Hz (required)
            %   params : Struct of parameter overrides (optional)
            %   tProbe : timeProbe instance (optional)

            if nargin < 1 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('hebbASR:invalidSrate', 'srate must be a positive numeric scalar.');
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

            [obj.B, obj.A] = hebbASR.buildSpectralFilter(srate);
            obj.modifiedMask = false(1, 0);
        end

        % ── Calibration ──────────────────────────────────────────────
        function calibrate(obj, X)
            % CALIBRATE - HMO-ASR baseline, then seed Hebbian weights.
            %   Faithful to calibrate_hebb_asr.m:
            %     1) Full HMO-ASR calibration (z_dev gate, Welford, early-eigen)
            %     2) Derive learning rates from eta and tau
            %     3) Seed W = E.' and M = I (rows of M\W span principal subspace)

            if isempty(X) || ~isnumeric(X)
                error('hebbASR:calibrate', 'Calibration data must be non-empty numeric.');
            end
            X = obj.ensureOrientation(X);

            mode = lower(strtrim(obj.hebb_mode));
            if ~ismember(mode, {'psp','psw'})
                error('hebbASR:calibrate', 'hebb_mode must be ''psp'' or ''psw''.');
            end
            obj.hebb_mode = mode;

            fprintf('--- [hebbASR/%s] Calibrating ---\n', mode);
            obj.tProbe.start('calibration');

            [C, S] = size(X);
            obj.nchans = C;
            X(~isfinite(X)) = 0;

            % Spectral shaping — preserve IIR state
            [Xf, iirstate_cal] = filter(obj.B, obj.A, double(X), [], 2);
            Xf = Xf.';   % [S x C]

            Nw  = max(1, round(obj.window_len * obj.srate));
            hop = max(1, round(Nw * (1 - obj.window_overlap)));
            win_starts = 1 : hop : max(1, S - Nw + 1);

            % HMO calibration accumulators
            n_tot = 0;  m_ch = zeros(C,1);  C_ch = zeros(C);
            E_prev = [];  sqrtC_run = eye(C);  mu_pc = [];  var_pc = [];
            n_dev = 0;  mu_dev = 0;  var_dev = 1;

            for ws = 1:numel(win_starts)
                idx = win_starts(ws) : win_starts(ws) + Nw - 1;
                if idx(end) > S, break; end

                W = Xf(idx, :);   % no Level-1 prune during calibration

                % Bootstrap first window
                if isempty(E_prev)
                    n_star = size(W,1);
                    if n_star < 2, continue; end
                    m_star = mean(W,1).';
                    Xc = W - m_star.';
                    C_ch = (Xc.'*Xc) / (n_star-1);
                    m_ch = m_star;  n_tot = n_star;
                    [E_prev, Dc] = eig((C_ch+C_ch.')/2 + 1e-10*eye(C));
                    lam = max(0, real(diag(Dc)));
                    sqrtC_run = E_prev * diag(sqrt(lam)) * E_prev.';
                    Y = (E_prev.' * W.').';
                    pc_rms = sqrt(mean(Y.^2,1)).';
                    mu_pc = pc_rms;  var_pc = ones(C,1)*1e-6;
                    continue;
                end

                % Whitened covariance deviation gate
                Cw    = cov(W, 1);
                Delta = sqrtC_run \ Cw / sqrtC_run;
                dev   = norm(Delta - eye(C), 'fro');

                if n_dev == 0
                    mu_dev = dev;  var_dev = 1e-6;  n_dev = 1;  z_dev = 0;
                else
                    z_dev = (dev - mu_dev) / sqrt(var_dev + 1e-12);
                end

                if z_dev < obj.z_clean_range(1) || z_dev > obj.z_clean_range(2)
                    continue;
                end
                [mu_dev, var_dev, n_dev] = hebbASR.updateMeanVar(mu_dev, var_dev, n_dev, dev);

                % Welford update
                n_star = size(W,1);
                if n_star < 2, continue; end
                m_star = mean(W,1).';
                Xc     = W - m_star.';
                C_star = (Xc.'*Xc) / (n_star-1);

                n_new = n_tot + n_star;
                delta = m_ch - m_star;
                m_new = m_ch + (m_star - m_ch)*(n_star/n_new);
                C_new = (((n_tot-1)*C_ch) + ((n_star-1)*C_star) + ...
                         (n_tot*n_star/n_new)*(delta*delta.')) / (n_new-1);
                n_tot = n_new;  m_ch = m_new;  C_ch = C_new;

                [E_prev, lam, mu_pc, var_pc] = hebbASR.refreshEigenbasis( ...
                    E_prev, C_ch, obj.eigen_threshold, mu_pc, var_pc);
                sqrtC_run = E_prev * diag(sqrt(lam)) * E_prev.';

                Y      = (E_prev.' * W.').';
                pc_rms = sqrt(mean(Y.^2,1)).';
                if isempty(mu_pc)
                    mu_pc = pc_rms;  var_pc = ones(C,1)*1e-6;
                end
                [mu_pc, var_pc, ~] = hebbASR.updateMeanVar(mu_pc, var_pc, ws, pc_rms, true);
            end

            % Fallbacks
            if isempty(E_prev) || n_tot == 0
                E_prev    = eye(C);
                lam       = max(1e-12, real(diag((Xf.'*Xf)/max(1,S))));
                sqrtC_run = E_prev * diag(sqrt(lam)) * E_prev.';
            end
            if isempty(mu_pc)
                lam_now = max(0, real(diag(E_prev.'*C_ch*E_prev)));
                mu_pc   = sqrt(lam_now);  var_pc = ones(C,1)*1e-6;
            end

            % Store HMO model
            obj.n_tot  = n_tot;   obj.m_ch  = m_ch;   obj.C_ch  = C_ch;
            obj.E      = E_prev;  obj.sqrtC = sqrtC_run;
            obj.mu_pc  = mu_pc;   obj.var_pc = var_pc;
            obj.Tpc    = mu_pc + obj.cutoff * sqrt(max(var_pc, 0));
            obj.n_rms  = 0;  obj.ch_mu = zeros(C,1);  obj.ch_var = ones(C,1);

            % Derive Hebbian learning rates
            switch mode
                case 'psp'
                    % PSP (Alg. 1): rW = 2η, rM = rW/(2τ)
                    obj.rW = 2 * obj.eta;
                    obj.rM = obj.rW / (2 * obj.tau);
                case 'psw'
                    % PSW (Alg. 2): rW = η, rM = rW/τ
                    obj.rW = obj.eta;
                    obj.rM = obj.rW / obj.tau;
            end

            % Seed Hebbian weights from HMO eigenbasis
            % Rows of M\W span the principal subspace → start consistent with E
            obj.W         = E_prev.';   % [C x C] feedforward
            obj.M         = eye(C);     % [C x C] lateral (identity = no whitening yet)
            obj.hebb_step = 0;

            % Preserve IIR state
            obj.iir_state = iirstate_cal;
            obj.carry     = [];

            obj.probe = hebbASRprobe(C);
            obj.probe.setBaseline(obj);

            obj.nsamples     = 0;
            obj.modifiedMask = false(1, 0);

            obj.tProbe.stop('calibration');
            fprintf('--- [hebbASR/%s] Calibration complete ---\n', mode);
        end

        % ── Process ──────────────────────────────────────────────────
        function cleanedSignal = process(obj, data)
            % PROCESS - Online artifact rejection with Hebbian adaptive model.
            %   Faithful to process_hebb_asr.m.

            if isempty(obj.E)
                error('hebbASR:notCalibrated', 'Run calibrate() before process().');
            end
            if isempty(data) || ~isnumeric(data)
                error('hebbASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            data = obj.ensureOrientation(data);

            obj.tProbe.start('process');
            [C, S] = size(data);

            mdims = obj.maxdims;
            if mdims > 0 && mdims < 1, mdims = round(C * mdims); end

            Nw  = max(1, round(obj.window_len * obj.srate));
            hop = max(1, round(obj.stepsize   * obj.srate));
            P   = max(0, round(obj.lookahead  * obj.srate));

            obj.modifiedMask(obj.nsamples + (1:S)) = false;

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
            Xf = Xf.';
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
                [Wpp, ~] = hebbASR.pcaEigZPrune(W_filt, 1.5);

                % Level-2 channel-RMS z-gate
                ch_rms = sqrt(mean(Wpp.^2, 1)).';
                [obj.ch_mu, obj.ch_var, obj.n_rms] = ...
                    hebbASR.updateMeanVar(obj.ch_mu, obj.ch_var, obj.n_rms, ch_rms);
                z_ch     = (ch_rms - obj.ch_mu) ./ max(1e-12, sqrt(obj.ch_var));
                is_clean = all(z_ch >= obj.z_clean_range(1) & z_ch <= obj.z_clean_range(2));
                obj.tProbe.stop('gating');

                % Clean-model update: Welford + Hebbian + derive E
                if is_clean
                    obj.tProbe.start('model_update');
                    obj.updateHebbianModel(Wpp, w);
                    obj.tProbe.stop('model_update');
                end

                % Reconstruction (identical projector logic to hmoASR)
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

            obj.carry     = data_ext(:, end-P+1:end);
            cleanedSignal = cleaned(:, P+1:end);
            obj.nsamples  = obj.nsamples + S;

            obj.tProbe.stop('process');
        end

    end % public methods

    % ================================================================
    methods (Access = private)

        function X = ensureOrientation(obj, X)
            if ~isnumeric(X) || isempty(X)
                error('hebbASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            if ~isempty(obj.nchans)
                if size(X,1) ~= obj.nchans && size(X,2) == obj.nchans
                    X = X.';
                elseif size(X,1) ~= obj.nchans
                    error('hebbASR:channelMismatch', ...
                        'Expected %d channels, got [%d x %d].', ...
                        obj.nchans, size(X,1), size(X,2));
                end
            else
                if size(X,1) > size(X,2), X = X.'; end
            end
        end

        function updateHebbianModel(obj, Wpp, win_idx)
            % UPDATEHEBBIANMODEL - Welford + PSP/PSW update + E from Hebbian weights.
            %   Faithful to process_hebb_asr.m sections (A)-(E).

            % (A) Welford covariance update
            n_star = size(Wpp,1);
            if n_star < 2, return; end

            m_star = mean(Wpp,1).';
            Xc     = Wpp - m_star.';
            C_star = (Xc.'*Xc) / (n_star-1);

            n_new = obj.n_tot + n_star;
            if obj.n_tot == 0
                m_new = m_star;  C_new = C_star;
            else
                delta = obj.m_ch - m_star;
                m_new = obj.m_ch + (m_star - obj.m_ch)*(n_star/n_new);
                C_new = (((obj.n_tot-1)*obj.C_ch) + ((n_star-1)*C_star) + ...
                         (obj.n_tot*n_star/n_new)*(delta*delta.')) / (n_new-1);
            end
            obj.n_tot = n_new;  obj.m_ch = m_new;  obj.C_ch = C_new;

            % (B) Hebbian / anti-Hebbian weight update
            Xk = Wpp.';   % [C x Nw]
            obj.hebb_step = obj.hebb_step + 1;

            switch obj.hebb_mode
                case 'psp'
                    [obj.W, obj.M] = hebbASR.hebbUpdatePSP(obj.W, obj.M, Xk, obj.rW, obj.rM);
                case 'psw'
                    [obj.W, obj.M] = hebbASR.hebbUpdatePSW(obj.W, obj.M, Xk, obj.rW, obj.rM);
            end

            % (C) Derive eigenbasis from Hebbian weights
            %     Rows of M\W span the principal subspace (Tsai et al. Eq. 6)
            Psub   = (obj.M \ obj.W).';   % [C x C] — columns span principal subspace
            [E_new, ~] = qr(Psub, 0);    % economy QR to orthonormalise
            obj.E  = E_new;

            % (D) Update sqrtC from global clean cov in Hebbian basis
            lam_c   = max(0, real(diag(obj.E.' * obj.C_ch * obj.E)));
            obj.sqrtC = obj.E * diag(sqrt(lam_c)) * obj.E.';

            % (E) Per-PC RMS stats & thresholds
            C = obj.nchans;
            Y      = (obj.E.' * Wpp.').';
            pc_rms = sqrt(mean(Y.^2,1)).';
            if isempty(obj.mu_pc)
                obj.mu_pc  = pc_rms;
                obj.var_pc = ones(C,1)*1e-6;
            end
            [obj.mu_pc, obj.var_pc, ~] = hebbASR.updateMeanVar( ...
                obj.mu_pc, obj.var_pc, win_idx, pc_rms, true);
            obj.Tpc = obj.mu_pc + obj.cutoff * sqrt(max(obj.var_pc, 0));
        end

        function [Rproj, Ek, lam, Tproj, reject] = buildProjector(obj, W_filt, Nw, mdims)
            COVw     = (W_filt.' * W_filt) / max(1, Nw);
            [Ek, Dk] = eig((COVw + COVw.')/2);
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
                C1s    = (C1+C1.')/2;
                [V, D] = eig(C1s);
                d      = diag(D);  d(d < eps) = eps;
                i_sqrt = V * diag(1./sqrt(d)) * V.';
                mat    = i_sqrt * C2 * i_sqrt;
                vals   = eig((mat+mat.')/2);
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
                [B, A] = butter(4, 2*[0.5 min(45, srate/2-2)]/srate, 'bandpass');
            end
        end

        function [W_new, M_new] = hebbUpdatePSP(W, M, Xk, rW, rM)
            % PSP — Algorithm 1, Tsai et al. (2024)
            Yk    = M \ (W * Xk);
            W_new = (1 - rW)*W + rW * (Yk * Xk.');
            M_new = (1 - rM)*M + rM * (Yk * Yk.');
        end

        function [W_new, M_new] = hebbUpdatePSW(W, M, Xk, rW, rM)
            % PSW — Algorithm 2, Tsai et al. (2024)
            C     = size(W, 1);
            Yk    = M \ (W * Xk);
            W_new = (1 - rW)*W + rW * (Yk * Xk.');
            M_new = (M - rM * eye(C)) + rM * (Yk * Yk.');
        end

        function [E_new, lam, mu_pc, var_pc] = refreshEigenbasis(E_prev, C_ch, thr, mu_pc, var_pc)
            if isempty(E_prev)
                need = true;
            else
                Rtest = E_prev.' * C_ch * E_prev;
                offd  = sum(abs(Rtest(:))) - sum(abs(diag(Rtest)));
                diagv = sum(abs(diag(Rtest))) + 1e-12;
                need  = (offd/diagv) >= thr;
            end
            if need
                E_old = E_prev;
                if isempty(E_old), E_old = eye(size(C_ch,1)); end
                [E_new, Dc] = eig((C_ch+C_ch.')/2);
                lam = max(0, real(diag(Dc)));
                if ~isempty(mu_pc)
                    [mu_pc, var_pc] = hebbASR.projectStatsNewBasis(mu_pc, var_pc, E_old, E_new);
                end
            else
                E_new = E_prev;
                lam   = max(0, real(diag(E_new.'*C_ch*E_new)));
            end
        end

        function [Wclean, kept] = pcaEigZPrune(W, zcut)
            COV  = (W.'*W) / max(1, size(W,1));
            [E, D] = eig((COV+COV.')/2);
            lam  = real(diag(D));
            z    = (lam - mean(lam)) / (std(lam)+1e-12);
            kept = z <= zcut;
            E2   = E(:, kept);
            if isempty(E2)
                [~,i] = min(lam);  E2 = E(:,i);
                kept = false(size(kept));  kept(i) = true;
            end
            Wclean = (E2 * (E2.'*W.')).';
        end

        function [mu, varv, n] = updateMeanVar(mu, varv, n, x, robust_init)
            if nargin < 5, robust_init = false; end
            x = x(:);
            if n == 0
                mu = x;  varv = ones(size(x))*(robust_init*1e-6);  n = 1;  return;
            end
            n_new  = n + 1;
            delta  = x - mu;
            mu_new = mu + delta./n_new;
            M2_old = varv.*max(1, n-1);
            M2_new = M2_old + delta.*(x - mu_new);
            varv   = max(M2_new./max(1, n_new-1), 0);
            mu = mu_new;  n = n_new;
        end

        function [mu_new, var_new] = projectStatsNewBasis(mu_old, var_old, E_old, E_new)
            Q       = E_old.' * E_new;
            s2_old  = mu_old.^2 + var_old;
            s2_new  = (Q.^2).' * s2_old;
            var_new = (Q.^2).' * var_old;
            mu_new  = sqrt(max(s2_new - var_new, 0));
        end

    end % static private methods

end
