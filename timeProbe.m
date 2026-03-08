classdef timeProbe < handle
    % timeProbe - Generic high-resolution performance profiling utility.
    %   Provides O(1) timing measurements via pre-allocated memory buffers.
    %   Designed for high-frequency loops where dynamic allocation overhead
    %   would bias execution timing.
    %
    %   Usage:
    %       tp = timeProbe();           % default 1e5 capacity per stage
    %       tp = timeProbe(5e4);        % custom capacity
    %
    %       tp.start('calibration');
    %       % ... work ...
    %       tp.stop('calibration');
    %       tp.stop('reconstruction', struct('trivial', true));
    %
    %       stats = tp.summarize('calibration');
    %       all   = tp.summarizeAll();
    %       tp.listStages();

    properties (Access = public)
        enabled  logical = true     % Global toggle — set false to zero profiling overhead
    end

    properties (SetAccess = private, GetAccess = public)
        % FIX #9: capacity is readable but only settable in the constructor
        capacity (1,1) double = 1e5
    end

    properties (Access = private)
        timers  struct = struct()   % Active tic identifiers keyed by stage name
        records struct = struct()   % Pre-allocated record arrays keyed by stage name
        indices struct = struct()   % Next-write index per stage
        warned  struct = struct()   % FIX #2: tracks whether overflow warning was issued
    end

    methods

        function obj = timeProbe(capacity)
            % Constructor
            %   capacity: (optional) max records per stage. Default 1e5.
            % FIX #4: removed redundant else-branch default assignment
            if nargin > 0
                obj.capacity = capacity;
            end
        end

        function start(obj, name)
            % START - Begin timing a named stage.
            if ~obj.enabled, return; end
            obj.timers.(name) = tic;
        end

        function stop(obj, name, meta)
            % STOP - Record elapsed time for a named stage. O(1).
            %   name : Stage name string — must match a prior start() call.
            %   meta : (optional) scalar struct of metadata for this record.
            if ~obj.enabled, return; end

            % FIX #1: guard against stop() called before start()
            if ~isfield(obj.timers, name) || isempty(obj.timers.(name))
                error('timeProbe:noTimer', ...
                    'stop() called for "%s" but start() was never called.', name);
            end

            dt = toc(obj.timers.(name));
            obj.timers.(name) = [];   % clear so a repeat stop() is caught

            % Lazy buffer initialisation on first call for this stage
            if ~isfield(obj.records, name)
                obj.initBuffer(name);
            end

            currIdx = obj.indices.(name);

            % FIX #2: warn once when buffer fills rather than silently dropping
            if currIdx > obj.capacity
                if ~isfield(obj.warned, name) || ~obj.warned.(name)
                    warning('timeProbe:bufferFull', ...
                        'Buffer full for stage "%s" (capacity=%d). ' ,...
                        'Further records are dropped. Increase capacity if needed.', ...
                        name, obj.capacity);
                    obj.warned.(name) = true;
                end
                return;
            end

            obj.records.(name)(currIdx).time = dt;
            if nargin > 2
                obj.records.(name)(currIdx).meta = meta;
            end
            obj.indices.(name) = currIdx + 1;
        end

        function stats = summarize(obj, name, metaField, metaValue)
            % SUMMARIZE - Return timing statistics for a named stage.
            %   summarize(name)                       — all records
            %   summarize(name, metaField, metaValue) — filtered by metadata
            %
            %   Returns struct with fields: mean, std, min, max, n (seconds)

            if ~isfield(obj.records, name)
                error('timeProbe:unknownStage', ...
                    'No records found for stage "%s". Known stages: %s', ...
                    name, strjoin(obj.listStages(), ', '));
            end

            n    = obj.indices.(name) - 1;
            data = obj.records.(name)(1:n);

            if nargin < 3
                stats = obj.calc([data.time]);
            else
                % FIX #5: note — isfield([], f) returns false so the && short-circuit
                % protects against records with empty meta. Behaviour is correct.
                mask = arrayfun(@(x) isfield(x.meta, metaField) && ...
                                     isequal(x.meta.(metaField), metaValue), data);

                % FIX #6: warn when filter matches nothing but records exist
                if n > 0 && sum(mask) == 0
                    warning('timeProbe:noMatch', ...
                        'summarize("%s"): %d records exist but none match %s==%s.', ...
                        name, n, metaField, mat2str(metaValue));
                end
                stats = obj.calc([data(mask).time]);
            end
        end

        function result = summarizeAll(obj)
            % SUMMARIZEALL - Return timing statistics for every recorded stage.
            %   Returns a struct where each field is a stage name containing
            %   the same stats struct as summarize().
            %
            % FIX #8: new method — avoids calling summarize() per stage in loops
            stages = obj.listStages();
            result = struct();
            for i = 1:numel(stages)
                result.(stages{i}) = obj.summarize(stages{i});
            end
        end

        function stages = listStages(obj)
            % LISTSTAGES - Return a cell array of all recorded stage names.
            % FIX #7: new method — essential for debugging multi-stage experiments
            stages = fieldnames(obj.records);
        end

        function reset(obj)
            % RESET - Clear all recorded data and reset all indices.
            % FIX #3: now clears records AND warned in addition to indices/timers,
            % so no stale data remains in the buffer tails after reset.
            obj.timers  = struct();
            obj.records = struct();
            obj.indices = struct();
            obj.warned  = struct();
        end

    end % end public methods

    methods (Access = private)

        function initBuffer(obj, name)
            % Pre-allocate a struct array of fixed capacity for this stage.
            % Template trick avoids O(N) growth during the experiment loop.
            template(obj.capacity) = struct('time', 0, 'meta', []);
            obj.records.(name)     = template;
            obj.indices.(name)     = 1;
            obj.warned.(name)      = false;
        end

        function s = calc(~, vals)
            % CALC - Compute summary statistics over a vector of elapsed times.
            if isempty(vals)
                s.mean = NaN;
                s.std  = NaN;
                s.min  = NaN;
                s.max  = NaN;
                s.n    = 0;
            else
                s.mean = mean(vals);
                s.std  = std(vals);
                s.min  = min(vals);
                s.max  = max(vals);
                s.n    = numel(vals);
            end
        end

    end % end private methods

end
