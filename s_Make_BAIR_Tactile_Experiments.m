% Make Tactile Experiment Files


% Prompt for ExperimentSpecs
[experimentSpecs, whichSite, selectionMade] = bairExperimentSpecs('prompt', true);
if ~selectionMade, return; end

% Which experiment to make?
[experimentType, selectionMade] = bairWhichExperimentTactile();
if ~selectionMade, return; end

% Set some defaults for all the experiments
TR              = 0.850;      % seconds
stimDiameterDeg = 16.6;       % degrees

% Generate stimulus template and set some defaults
if whichSite == 3
    experimentSpecs(3,1) = {'BAIR_ACER'};
end

stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);
stimParams.NIdaqRate = 1000;
stimParams.NIdaqNames = 'cDAQ1mod1';

% make a figure showing the stimulus pattern
makeFigure = 'true'; 

% Find the selected experiment
switch experimentType

    case 'SIMPLEHANDSWEEP' %simple hand stimulation in different patterns for now
        
        directions = {'ascending', 'descending', 'random', 'staircase1','staircase2'}; %makes one of each by default
        stimDurationSeconds = 6; % seconds
        numberOfStimulators = 5;
        numberOfRuns = 1; 
        
        for runNum = 1:numberOfRuns
            % MAKE TASK EXPERIMENT
            stimMakeTactileHandExperiment(stimParams,  runNum, stimDurationSeconds,...
                TR, numberOfStimulators,directions ,makeFigure)
        end 
        
    case 'BlOCKED'
        
        numStimulatorsInBlock = 2;
        directions = {'ascending', 'descending', 'random', 'staircase1','staircase2'}; %makes one of each by default
        stimDurationSeconds = 6; % seconds
        numberOfStimulators = 10;
        numberOfRuns = 1; 
        
        %stimMakeTactileBlockedExperiment(stimParams,  runNum, stimDurationSeconds, ...
        % TR, numberOfStimulators,numStimulatorsInBlock,directions ,makeFigure)
end



