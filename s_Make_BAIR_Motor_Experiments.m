% Make Motor Experiment Files


% Prompt for ExperimentSpecs
[experimentSpecs, whichSite, selectionMade] = bairExperimentSpecs('prompt', true);
if ~selectionMade, return; end

% Which experiment to make?
[experimentType, selectionMade] = bairWhichExperiment();
if ~selectionMade, return; end

% Set some defaults for all the experiments
TR              = 0.850;      % ms
stimDiameterDeg = 16.6;       % degrees

% Generate stimulus template and set some defaults
stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);

% Find the selected experiment
switch experimentType

    
     case 'BOLDHAND'
        stimDurationSeconds = 0.500; % seconds
        onsetTimeMultiple   = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR
        numberOfRuns = 1; 
        
        for runNum = 1:numberOfRuns
            % MAKE TASK EXPERIMENT
            stimMakeBoldHandExperiment(stimParams,  runNum, stimDurationSeconds, onsetTimeMultiple, TR)
        end 
        
        
    case {'GESTURES' 'GESTURESTRAINING'}
        stimDurationSeconds    = 5;
        numberOfRuns = 1;
        
        for runNum = 1:numberOfRuns
            % MAKE TASK EXPERIMENT
            stimMakeGesturesExperiment(stimParams, runNum, TR, stimDurationSeconds);
        end
        
    case {'FINGERMAPPINGLEFT', 'FINGERMAPPINGRIGHT'}
        numberOfRuns = 1;
        for runNum = 1:numberOfRuns
            stimMakeFingermappingExperiment(stimParams,  runNum, TR, experimentType)
        end
               
    case 'BOLDSAT'
        
end

