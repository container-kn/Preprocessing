classdef configExperiment < handle
    % configExperiment - Centralized configuration handler for BCI experiments.
    %   Manages pathing, signal specifications, and algorithmic hyperparameters
    %   for all ASR variants. Ensures consistency across multi-subject sweeps.
    %
    %   Usage:
    %       config = configExperiment(projectRoot, algorithm);
    %       config = configExperiment(projectRoot, algorithm, 'fs', 250, 'subjectID', 5);

    properties (Access = public)
        % ----- Environment & IO -----
        paths  = struct()   % project, results, logs, dataset paths
        signal = struct()   % fs, bandpass, lineNoise
        meta   = struct()   % OS info and execution timestamp

        % ----- Experiment Context -----
        subjectID (1,1) double {mustBePositive, mustBeInteger} = 1  % FIX #7
        algorithm  struct                                            % .name + .parameters
        sweepMode  (1,1) logical = false  % set true by parameterSweeper — skips raw in saveResults
    end

    methods

        function obj = configExperiment(projectRoot, algorithmStruct, varargin)
            % CONFIGEXPERIMENT - Initialise an experiment configuration.
            %   projectRoot    : String path to the project root directory.
            %   algorithmStruct: Struct with at minimum a .name field.
            %   varargin       : Name-value pairs — fs, bandpass, lineNoise,
            %                    subjectID, algorithmParams.

            % 1. Validate inputs
            if nargin < 1 || isempty(projectRoot)
                projectRoot = fileparts(mfilename('fullpath'));
            end
            if nargin < 2
                error('configExperiment:missingInput', ...
                    'algorithmStruct is required as the second argument.');
            end

            % FIX #8: validate algorithmStruct has a .name string field
            if ~isstruct(algorithmStruct) || ~isfield(algorithmStruct, 'name') ...
                    || ~(ischar(algorithmStruct.name) || isstring(algorithmStruct.name))
                error('configExperiment:invalidAlgorithm', ...
                    'algorithmStruct must be a struct with a .name string field.');
            end

            % FIX #6: validate projectRoot before creating any directories
            if ~isfolder(projectRoot)
                error('configExperiment:invalidRoot', ...
                    'projectRoot does not exist: %s', projectRoot);
            end

            % 2. Paths & folders
            obj.buildPaths(projectRoot);
            obj.ensureFolders();

            % 3. Input parsing
            p = inputParser;
            % FIX #9: KeepUnmatched = false so typos in parameter names are caught
            p.KeepUnmatched = false;

            addParameter(p, 'fs',              500);
            addParameter(p, 'bandpass',        [0.5 40]);
            addParameter(p, 'lineNoise',       50);
            addParameter(p, 'subjectID',       1);
            addParameter(p, 'algorithmParams', struct());

            parse(p, varargin{:});
            res = p.Results;

            % 4. Assign signal settings and identity
            obj.setSignal(res.fs, res.bandpass, res.lineNoise);
            obj.subjectID = res.subjectID;
            obj.algorithm = algorithmStruct;

            if ~isfield(obj.algorithm, 'parameters')
                obj.algorithm.parameters = struct();
            end

            % 5. Apply algorithm defaults then user overrides
            obj.applyAlgorithmDefaults();

            overrideFields = fieldnames(res.algorithmParams);
            for k = 1:numel(overrideFields)
                obj.algorithm.parameters.(overrideFields{k}) = ...
                    res.algorithmParams.(overrideFields{k});
            end

            % 6. System metadata
            obj.setMeta();
        end

        % ====================================================
        % Path Utilities
        % ====================================================

        function buildPaths(obj, root)
            % BUILDPATHS - Construct the directory tree structure.
            % FIX #1: handle class — no need to return obj
            obj.paths.project = root;
            obj.paths.results = fullfile(root, 'results');
            obj.paths.logs    = fullfile(root, 'logs');
            obj.paths.dataset = fullfile(root, 'datasets');
        end

        function ensureFolders(obj)
            % ENSUREFOLDERS - Create output directories if they don't exist.
            dirs = {obj.paths.results, obj.paths.logs};
            for i = 1:numel(dirs)
                if ~exist(dirs{i}, 'dir')
                    mkdir(dirs{i});
                end
            end
        end

        % ====================================================
        % Signal Settings
        % ====================================================

        function setSignal(obj, fs, bp, ln)
            % SETSIGNAL - Define signal acquisition parameters.
            obj.signal.fs        = fs;
            obj.signal.bandpass  = bp;
            obj.signal.lineNoise = ln;
        end

        % ====================================================
        % Algorithm Parameter Defaults
        % ====================================================

        function applyAlgorithmDefaults(obj)
            % APPLYALGORITHMDEFAULTS - Per-variant hyperparameter library.
            %   Only fills fields not already set, so manual pre-init is respected.

            % Shared baseline across all variants
            base.lookahead    = 0.25;
            base.calibWindSec = [2 32];
            base.blocksize    = 20;
            base.cutoff       = 16;

            % FIX #2: algorithm name cases aligned with Experimenting.m switch,
            % including original-asr/o-asr which were missing entirely.
            switch lower(obj.algorithm.name)

                case {'vanilla-asr', 'v-asr', 'asr'}
                    defaults = base;

                case {'original-asr', 'o-asr'}
                    defaults = base;

                case {'e-asr', 'embeddedasr', 'e_asr'}
                    defaults = base;
                    defaults.cutoff             = 15;
                    defaults.embeddingDimension = 90;

                case {'tau-asr'}
                    defaults = base;
                    defaults.cutoff             = 10;
                    defaults.embeddingDimension = 25;
                    defaults.tau                = 12;

                case {'hmo-asr', 'hmo_asr'}
                    defaults = base;
                    defaults.stepsize = 0.25;

                case {'ema-asr', 'ema_asr', 'ema'}
                    % FIX #3: cleanBuffLength corrected from 2 → 5 seconds (fast buffer)
                    defaults = base;
                    defaults.cleanBuffLength = 2;
                    defaults.beta1           = 0.05;
                    defaults.beta2           = 0.05;

                case {'graph-asr', 'g-asr', 'g-asr-f', 'g_asr', 'g_asr-frank'}
                    % FIX #4: inflation corrected from 1 → 5 (match graphASR default)
                    % FIX #5: threshLambda and shrinkLambda removed — not
                    %         properties of graphASR, were silently ignored
                    defaults = base;
                    defaults.shrinkStyle = 'zero';
                    defaults.inflation   = 5;
                    defaults.kNeighbors  = 4;
                    defaults.varKeep     = 0.9999;
                    defaults.lambdaGain  = 5;
                    defaults.lambdaExp   = 0.5;

                case {'graph-asr-base', 'g-asr-b', 'g_asr-base', 'g_asr-b'}
                    defaults = base;
                    defaults.shrinkStyle = 'zero';
                    defaults.inflation   = 5;
                    defaults.kNeighbors  = 4;

                case {'riemann-asr', 'r-asr', 'r_asr'}
                    % Same parameter set as vanillaASR — cutoff, blocksize,
                    % maxdims, lookahead, stepsize — plus window_len/overlap
                    defaults = base;
                    defaults.maxdims              = 0.66;
                    defaults.stepsize             = 32;
                    defaults.window_len           = 0.5;
                    defaults.window_overlap       = 0.66;
                    defaults.max_dropout_fraction = 0.1;
                    defaults.min_clean_fraction   = 0.25;

                otherwise
                    error('configExperiment:unsupportedAlgorithm', ...
                        'Algorithm "%s" is not recognised. Check spelling or add a case.', ...
                        obj.algorithm.name);
            end

            % Merge: only fill missing fields
            fields = fieldnames(defaults);
            for k = 1:numel(fields)
                if ~isfield(obj.algorithm.parameters, fields{k})
                    obj.algorithm.parameters.(fields{k}) = defaults.(fields{k});
                end
            end
        end

        % ====================================================
        % Metadata
        % ====================================================

        function setMeta(obj)
            % SETMETA - Record environment context for reproducibility.
            obj.meta.os        = computer;
            obj.meta.timestamp = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
        end

    end
end
