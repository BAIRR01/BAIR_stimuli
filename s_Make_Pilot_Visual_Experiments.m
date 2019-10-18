% Make Pilot Visual Experiment Files
%
% 1. Six category localizer
% 2. Eight category localizer
% 3. Kalanit category localizer
% 4. Scene face lateral
% 4. Kravitz scene categories
% 5. Bonner scene affordance

% Prompt for ExperimentSpecs
[experimentSpecs, whichSite, selectionMade] = pilotExperimentSpecs('prompt', true);
if ~selectionMade, return; end

% Which experiment to make?
[experimentType, selectionMade] = pilotWhichExperimentVisual();
if ~selectionMade, return; end

% Generate stimulus template

%   max stimulus radius (in deg)
%       16.6º is the height of the screen for the 3T display at Utrecht,
%       which is the smallest FOV among NYU-3T, UMC-3T, NYU-ECoG, UMC-ECoG,
%       NYU-MEG

TR              = 0.850;      % ms
stimDiameterDeg = 16.6;       % degrees

stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);

switch experimentType
    case {'SIXCATLOC'} %
        % Make SIXCATLOC experiment

        % We have 2 unique Master runs with fixed stimulus orders that will
        % be identical for other modalities. Timing (ISI) is different
        % between modalities; the fixation sequence is generated anew for
        % each new experiment (run)
        numberOfRuns           = 2;        
        onsetTimeMultiple      = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR (fMRI experiments only)
        
        for runNum = 1:numberOfRuns
            stimMakeLocalizerExperiment(stimParams, runNum, experimentType, onsetTimeMultiple, TR);
        end   
        
	case {'OBJECTDETECTION'} %
        % Make SIXCATLOC experiment

        % We have 2 unique Master runs with fixed stimulus orders that will
        % be identical for other modalities. Timing (ISI) is different
        % between modalities; the fixation sequence is generated anew for
        % each new experiment (run)
        numberOfRuns           = 2;        
        onsetTimeMultiple      = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR (fMRI experiments only)
        
        for runNum = 1:numberOfRuns
            stimMakeSceneExperiment(stimParams, runNum, experimentType, onsetTimeMultiple, TR);
        end   
end

