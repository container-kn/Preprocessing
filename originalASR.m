classdef originalASR < handle
    % originalASR - Wrapper for standard EEGLAB/BCILAB ASR functions.
    %   This class encapsulates asr_calibrate and asr_process while maintaining
    %   consistent object-oriented state management compatible with the other
    %   ASR variant classes (vanillaASR, emaASR, graphASR).
    %
    %   Constructor signature matches the other variants:
    %       obj = originalASR(X, srate, cutoff_val, timeProbe)

    properties (Access = public)
        % Core Dimensions
        nchans
        nsamples = 0        % FIX #4: tracks processed samples only, not calibration data
        modifiedMask

        % Algorithmic Parameters
        srate = 500
        cutoff = 20
        blocksize = 10      % FIX #7: now wired into asr_calibrate
        window_len = 0.5
        window_overlap = 0.66
        max_dropout_fraction = 0.1
        min_clean_fraction = 0.25
        maxdims = 0.66
        lookahead = 0.25    % NOTE: asr_process may interpret this as seconds or samples
                            % depending on BCILAB version — verify against your build.
        stepsize = 32

        % Calibration State (Standard ASR structure from asr_calibrate)
        state

        % Diagnostic Probes
        timeProbe = []
    end

    % FIX #5, #6: removed dead iir_state and carry properties — asr_process
    % manages both internally via the state struct.

    methods

        function obj = originalASR(X, srate, cutoff_val, timeProbe)
            % Constructor
            %   X          : Data matrix used to determine channel count [chans x samples]
            %                (calibration data — not stored, just sized)
            %   srate      : Sampling rate in Hz
            %   cutoff_val : SD threshold for artifact detection (default 20)
            %   timeProbe  : Timing probe object (required by harness)

            % FIX #10: Input validation
            if nargin < 2 || isempty(srate) || ~isnumeric(srate) || srate <= 0
                error('originalASR: srate must be a positive numeric scalar.');
            end
            if isempty(X) || ~isnumeric(X)
                error('originalASR: X must be a non-empty numeric matrix.');
            end

            obj.srate = srate;

            if nargin >= 3 && ~isempty(cutoff_val)
                obj.cutoff = cutoff_val;
            end

            % FIX #2: timeProbe is now a proper constructor argument
            if nargin >= 4 && ~isempty(timeProbe)
                obj.timeProbe = timeProbe;
            end

            X            = obj.determineSignalShape(X);
            obj.nchans   = size(X, 1);
            % FIX #4: nsamples starts at 0 — it tracks processed samples, not
            % calibration data size. The two are conceptually separate.
            obj.nsamples    = 0;
            obj.modifiedMask = false(1, 0);
        end

        function X = determineSignalShape(obj, X)
            % Orient data to [chans x samples].
            % FIX #9: post-calibration, use nchans to enforce orientation rather
            % than guessing — avoids silent errors on square matrices.
            if ~isnumeric(X) || isempty(X)
                error('originalASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end
            if ~isempty(obj.nchans)
                if size(X, 1) == obj.nchans
                    return;
                elseif size(X, 2) == obj.nchans
                    X = X';
                else
                    error('originalASR:channelMismatch', ...
                        'Data has %d rows and %d cols; expected %d channels.', ...
                        size(X,1), size(X,2), obj.nchans);
                end
            else
                % Pre-calibration heuristic
                if size(X, 1) > size(X, 2)
                    X = X';
                end
            end
        end

        function calibrate(obj, X)
            % calibrate - Establish baseline via asr_calibrate().
            %
            % FIX #10: input validation
            if isempty(X) || ~isnumeric(X)
                error('originalASR:calibrate', 'Calibration data must be a non-empty numeric matrix.');
            end

            X = obj.determineSignalShape(X);

            obj.timeProbe.start('calibration');
            % FIX #7: blocksize is now passed as the 4th argument to asr_calibrate.
            % Argument order: data, srate, cutoff, blocksize, window_len,
            %                 window_overlap, max_dropout_fraction, min_clean_fraction
            obj.state = asr_calibrate(double(X),        ...
                obj.srate,                              ...
                obj.cutoff,                             ...
                obj.blocksize,                          ...
                obj.window_len,                         ...
                obj.window_overlap,                     ...
                obj.max_dropout_fraction,               ...
                obj.min_clean_fraction);
            obj.timeProbe.stop('calibration');

            % Reset runtime counters on re-calibration
            obj.nsamples     = 0;
            obj.modifiedMask = false(1, 0);
        end

        function cleanedSignal = process(obj, data)
            % process - Clean data via asr_process().

            % FIX #10: guard against uncalibrated use
            if isempty(obj.state)
                error('originalASR:notCalibrated', 'Run calibrate() before process().');
            end
            if isempty(data) || ~isnumeric(data)
                error('originalASR:invalidInput', 'Data must be a non-empty numeric matrix.');
            end

            data = obj.determineSignalShape(data);

            obj.timeProbe.start('process');
            [cleanedSignal, obj.state] = asr_process(double(data), ...
                obj.srate,          ...
                obj.state,          ...
                obj.window_len,     ...
                obj.window_overlap, ...
                obj.maxdims,        ...
                obj.stepsize,       ...
                obj.lookahead);
            obj.timeProbe.stop('process');  % FIX #1: was .start() — timer never closed

            % FIX #3: derive modifiedMask by comparing input vs output.
            % Samples where the signal was altered = reconstructed by ASR.
            % Tolerance of 1e-10 accounts for floating-point pass-through noise.
            S          = size(data, 2);
            S_out      = size(cleanedSignal, 2);
            % asr_process may return fewer samples than input (lookahead trim)
            n_compare  = min(S, S_out);
            diff_norm  = sum((data(:, 1:n_compare) - cleanedSignal(:, 1:n_compare)).^2, 1);
            new_mask   = diff_norm > 1e-10;

            % Extend modifiedMask to cover new samples
            obj.modifiedMask(obj.nsamples + (1:S_out)) = false;
            obj.modifiedMask(obj.nsamples + (1:n_compare)) = new_mask;

            obj.nsamples = obj.nsamples + S_out;
        end

    end
end
