% Make Pilot Visual Experiment Files
%
% 1. Six category localizer
% 2. Six category localizer temporal
% 3. Object detection
% Other options (To do):
% 4. Eight category localizer
% 5. Kalanit category localizer
% 6. Scene face lateral
% 7. Kravitz scene categories
% 8. Bonner scene affordance

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
    case {'SIXCATLOC', 'SIXCATLOCTEMPORAL', 'SIXCATLOCTEMPORALDIFF'} %
        % Make SIXCATLOC experiment

        % We have two unique sets of stimuli for even and odd runs.
        % Stimulus order is randomized across runs.
        % For temporal, assignment of temporal condition to image is fixed
        % with a seed based on runnumber.
        % Fixation sequence is generated anew for each run.
        numberOfRuns           = 6;        
        onsetTimeMultiple      = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR (fMRI experiments only)
        
        for runNum = 1:numberOfRuns
            stimMakeLocalizerExperiment(stimParams, runNum, experimentType, onsetTimeMultiple, TR);
        end   
     
    case {'SIXCATLOCISIDIFF'} %
        
        % We have two unique sets of stimuli for even and odd runs.
        % Stimulus order is randomized across runs.
        % For temporal, assignment of temporal condition to image is fixed
        % with a seed based on runnumber.
        % Fixation sequence is generated anew for each run.
        numberOfRuns           = 6;        
        onsetTimeMultiple      = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR (fMRI experiments only)
        
        for runNum = 1:numberOfRuns
            stimMakeLocalizerISIExperiment(stimParams, runNum, experimentType, onsetTimeMultiple, TR);
        end 
        
	case {'OBJECTDETECTION'} %
        % Make OBJECTDETECTION experiment

        % Stimulus order is randomized across runs.
        % Currently all 480 stimuli are loaded, subselection may be better.
        % Fixation sequence is generated anew for each run.
        numberOfRuns           = 2;        
        onsetTimeMultiple      = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR (fMRI experiments only)
        
        for runNum = 1:numberOfRuns
            stimMakeSceneExperiment(stimParams, runNum, experimentType, onsetTimeMultiple, TR);
        end   
        
	case {'SCENEFACELATERAL'} %
        % Make SCENEFACELATERAL experiment
        
        stimDiameterDeg = 20;       % degrees
        stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);

        % Stimulus order is randomized across runs.
        % Currently all 480 stimuli are loaded, subselection may be better.
        % Fixation sequence is generated anew for each run.
        numberOfRuns           = 2;        
        onsetTimeMultiple      = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR (fMRI experiments only)
        
        for runNum = 1:numberOfRuns
            stimMakeLateralVisualExperiment(stimParams, runNum, experimentType, onsetTimeMultiple, TR);
        end 
end

