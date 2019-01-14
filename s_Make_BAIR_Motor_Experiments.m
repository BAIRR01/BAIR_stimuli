% Make Motor Experiment Files


% Prompt for ExperimentSpecs
[experimentSpecs, whichSite, selectionMade] = bairExperimentSpecs('prompt', true);
if ~selectionMade, return; end

% Which experiment to make?
[experimentType, selectionMade] = bairWhichExperiment();
if ~selectionMade, return; end

% Set some defaults for all the experiments
TR              = 0.850;      % seconds
stimDiameterDeg = 16.6;       % degrees

% Generate stimulus template and set some defaults
stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);

% Find the selected experiment
switch experimentType

    case 'BOLDHAND'
        stimDurationSeconds = 0.500; % seconds
        %onsetTimeMultiple   = 0.170; % make the onsets multiple of 170 ms, which is 1/5 of the TR
        numberOfRuns = 1; 
        
        for runNum = 1:numberOfRuns
            % MAKE TASK EXPERIMENT
            stimMakeBoldHandExperiment(stimParams,  runNum, stimDurationSeconds, TR)
        end 
        
    case {'GESTURES' 'GESTURESPRACTICE' 'GESTURESLEARNING'}
        stimDurationSeconds    = 5; % seconds
        numberOfRuns = 1;
        
        for runNum = 1:numberOfRuns
            % MAKE TASK EXPERIMENT
            stimMakeGesturesExperiment(stimParams, runNum, TR, stimDurationSeconds, experimentType);
        end
        
    case {'FINGERMAPPINGLEFT', 'FINGERMAPPINGRIGHT'}
        numberOfRuns = 1;
        for runNum = 1:numberOfRuns
            stimMakeFingermappingExperiment(stimParams,  runNum, experimentType)
        end
               
    case 'BOLDSAT'
        stimDurationSeconds = 0.4; % seconds
        numberOfRuns = 1;
        movementRates = [1/3, .8, 1.3, 1.8]; % Hz 
        
        for runNum = 1:numberOfRuns
            for ii = 1:length(movementRates)
                movementRate = movementRates(ii);
                stimMakeBoldSatExperiment(stimParams, runNum, stimDurationSeconds, TR, movementRate,ii)
            end
        end
end

