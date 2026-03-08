classdef parameterSweeper < handle
    % parameterSweeper - Factorial parameter sweep engine for ASR benchmarking.
    %   Constructs a full-factorial grid from a parameter struct and runs one
    %   Experimenting pipeline per (subject, combination) pair, either serially
    %   or in parallel across combinations.
    %
    %   Usage:
    %       sweeper = parameterSweeper(baseConfig, algorithmStruct, [1 2 3 4]);
    %       sweeper.setSweep(struct('cutoff', [10 15 20], 'blocksize', [10 20]));
    %       sweeper.sweeprun();                    % serial
    %       sweeper.sweeprun('parallel', true);    % parallel over combinations

    properties (Access = public)
        baseConfig       % configExperiment instance (provides paths, signal settings)
        algorithmStruct  % Struct with .name field passed to configExperiment
        subjects         % Numeric array of subject IDs to sweep over
        paramGrid        % Struct of param → value vectors for factorial expansion
    end

    methods

        % ============================================
        % Constructor
        % ============================================
        function obj = parameterSweeper(baseConfig, algorithmStruct, subjects)
            if nargin < 3
                error('parameterSweeper:missingInput', ...
                    'Usage: parameterSweeper(baseConfig, algorithmStruct, subjects)');
            end
            obj.baseConfig      = baseConfig;
            obj.algorithmStruct = algorithmStruct;
            obj.subjects        = subjects(:)';   % enforce row vector
            obj.paramGrid       = struct();
        end

        % ============================================
        % Set sweep grid
        % ============================================
        function setSweep(obj, gridStruct)
            % SETSWEEP - Define the parameter grid for factorial expansion.
            %   gridStruct: struct where each field is a param name and value
            %               is a vector of values to sweep.
            %   Example:
            %       sweeper.setSweep(struct('cutoff', [10 15 20], 'blocksize', [10 20]));
            if ~isstruct(gridStruct) || isempty(fieldnames(gridStruct))
                error('parameterSweeper:invalidGrid', ...
                    'gridStruct must be a non-empty struct of param → value vectors.');
            end
            obj.paramGrid = gridStruct;
        end

        % ============================================
        % Main Sweep Runner
        % ============================================
        function sweeprun(obj, varargin)
            % SWEEPRUN - Execute the full factorial sweep.
            %   Options (name-value):
            %     'parallel' : logical, default false.
            %                  Parallelises over parameter combinations.

            % FIX #7: guard against empty grid
            if isempty(fieldnames(obj.paramGrid))
                error('parameterSweeper:emptyGrid', ...
                    'Call setSweep() before sweeprun().');
            end

            p = inputParser;
            p.KeepUnmatched = false;
            addParameter(p, 'parallel', false, @(x) islogical(x) || isnumeric(x));
            parse(p, varargin{:});
            useParallel = logical(p.Results.parallel);

            paramNames  = fieldnames(obj.paramGrid);
            paramValues = struct2cell(obj.paramGrid);

            % Unwrap: struct2cell may produce a cell-of-cells if the grid
            % was defined with {[...]} syntax (e.g. struct('cutoff',{[10 15]})).
            % Flatten each entry to a plain numeric vector before ndgrid.
            for i = 1:numel(paramValues)
                if iscell(paramValues{i})
                    paramValues{i} = [paramValues{i}{:}];
                end
            end

            % Build full factorial grid
            grids = cell(1, numel(paramValues));
            [grids{:}] = ndgrid(paramValues{:});
            totalComb  = numel(grids{1});
            nSubjects  = numel(obj.subjects);

            fprintf('\n=== Starting Sweep: %s ===\n', obj.algorithmStruct.name);
            fprintf('  Subjects            : %d\n', nSubjects);
            fprintf('  Parameter combos    : %d\n', totalComb);
            fprintf('  Total runs          : %d\n', nSubjects * totalComb);
            fprintf('  Parallel (combos)   : %d\n\n', useParallel);

            % FIX #6: parfor over combinations (better load balance),
            % serial over subjects inside — each combo is independent.
            % FIX #10: extract everything parfor needs as local variables
            % so the handle object itself isn't broadcast to workers.
            subjects_local      = obj.subjects;
            baseConfig_local    = obj.baseConfig;
            algStruct_local     = obj.algorithmStruct;

            if useParallel

                parfor k = 1:totalComb
                    override = parameterSweeper.buildOverride(paramNames, grids, k);

                    for s = 1:nSubjects
                        subjID = subjects_local(s);
                        % FIX #8: progress printed per worker (shows in parallel pool log)
                        fprintf('  [combo %d/%d] subject %d\n', k, totalComb, subjID);
                        parameterSweeper.runSingle( ...
                            baseConfig_local, algStruct_local, subjID, override);
                    end
                end

            else

                for k = 1:totalComb
                    override = parameterSweeper.buildOverride(paramNames, grids, k);

                    for s = 1:nSubjects
                        subjID = subjects_local(s);
                        % FIX #8: progress reporting in serial mode
                        runIdx = (k-1)*nSubjects + s;
                        fprintf('  [%d/%d] combo %d/%d — subject %d\n', ...
                            runIdx, nSubjects*totalComb, k, totalComb, subjID);
                        parameterSweeper.runSingle( ...
                            baseConfig_local, algStruct_local, subjID, override);
                    end
                end

            end

            fprintf('\n=== Sweep Completed ===\n');
        end

        % ============================================
        % Load results for a specific combo + subject
        % ============================================
        function results = loadResult(obj, subjID, override)
            % LOADRESULT - Load a saved result from disk for a given subject
            % and parameter override. Avoids re-running to inspect outputs.
            %
            % FIX #9: replaces the re-run pattern in the original nargout block.

            cfg = configExperiment( ...
                obj.baseConfig.paths.project, ...
                obj.algorithmStruct, ...
                'fs',              obj.baseConfig.signal.fs, ...
                'bandpass',        obj.baseConfig.signal.bandpass, ...
                'lineNoise',       obj.baseConfig.signal.lineNoise, ...
                'subjectID',       subjID, ...
                'algorithmParams', override);

            % Reconstruct the filename using the same pattern as saveResults()
            params   = cfg.algorithm.parameters;
            paramStr = parameterSweeper.buildParamString(params);
            fname    = sprintf('S%d_%s%s.mat', subjID, cfg.algorithm.name, paramStr);
            fpath    = fullfile(cfg.paths.results, fname);

            if ~isfile(fpath)
                error('parameterSweeper:fileNotFound', ...
                    'Result file not found: %s', fpath);
            end
            loaded  = load(fpath, 'results');
            results = loaded.results;
        end

    end % end public methods

    methods (Static, Access = private)

        % ============================================
        % Build override struct for combination k
        % ============================================
        function override = buildOverride(paramNames, grids, k)
            override = struct();
            for i = 1:numel(paramNames)
                override.(paramNames{i}) = grids{i}(k);
            end
        end

        % ============================================
        % Run a single (subject, override) experiment
        % ============================================
        function runSingle(baseConfig, algStruct, subjID, override)
            % RUNSINGLE - Build config, run experiment, save results.
            %   Static so it can be called cleanly from parfor without
            %   broadcasting the full parameterSweeper handle object.

            % FIX #3: configExperiment (correct class name)
            cfg = configExperiment( ...
                baseConfig.paths.project, ...
                algStruct, ...
                'fs',              baseConfig.signal.fs, ...
                'bandpass',        baseConfig.signal.bandpass, ...
                'lineNoise',       baseConfig.signal.lineNoise, ...
                'subjectID',       subjID, ...
                'algorithmParams', override);

            % Tell Experimenting we are in a sweep — raw signals saved once
            % separately via saveRawOnce(), not duplicated into every run file.
            cfg.sweepMode = true;

            % FIX #2: removed dead `probe = asrProbe()` line
            % FIX #4: Experimenting (correct class name)
            exp = Experimenting(cfg);

            % FIX #5: run() already calls saveResults() internally — don't call again
            exp.run();
        end

        function str = buildParamString(params)
            % Mirror of Experimenting.buildParamString for filename reconstruction.
            fields = fieldnames(params);
            str    = "";
            for i = 1:numel(fields)
                val = params.(fields{i});
                if iscell(val), val = val{1}; end
                if isnumeric(val)
                    if isscalar(val)
                        valStr = num2str(val);
                    else
                        valStr = sprintf('%g-', val);
                        valStr(end) = [];
                    end
                else
                    valStr = char(string(val));
                end
                str = str + "_" + fields{i} + "_" + valStr;
            end
        end

    end % end static private methods

end
