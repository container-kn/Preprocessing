function [signal,chanlabel,latency,chanlocs] = loadSubjectData(subjectNum, dataPath)
% loadSubjectData  Load EEG data for a given subject.
%
% Usage:
%   signal = loadSubjectData(1, 'H:\DynamicEASR\CodeBase\data\');
%
% Inputs:
%   subjectNum : integer, subject ID (e.g., 1 -> 'subject_1.mat')
%   dataPath   : char/string, path to the data folder
%
% Output:
%   signal     : [channels x samples] EEG data matrix

    % Build file path
    
    fileName = fullfile(dataPath, sprintf('subject_%d.mat', subjectNum));
    

    % Safety check
    if ~isfile(fileName)
        error('File not found: %s', fileName);
    end

    % Load .mat file
    dataStruct = load(fileName);

    % Expecting field "EEG.data"
    if ~isfield(dataStruct, 'EEG') || ~isfield(dataStruct.EEG, 'data')
        error('File %s does not contain EEG.data field.', fileName);
    end

    % Extract signal
    signal = dataStruct.EEG.data;
    chanlabel = {dataStruct.EEG.chanlocs.labels}; 
    latency = 150500;
    chanlocs = load_chanlocs_from_csv( ...
    fullfile(dataPath,'\headlocs.csv'));
    % chanlocs = importfile("H:\ArtifactSubspaceReconstruction\datasets\headlocs.csv");
end
