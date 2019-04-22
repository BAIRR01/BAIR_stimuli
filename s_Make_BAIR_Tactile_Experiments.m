% Make Tactile Experiment Files


% Prompt for ExperimentSpecs
[experimentSpecs, whichSite, selectionMade] = bairExperimentSpecs('prompt', true);
if ~selectionMade, return; end

% Which experiment to make?
[experimentType, selectionMade] = bairWhichExperimentTactile();
if ~selectionMade, return; end

% Set some defaults for all the experiments
TR              = 0.85;   % seconds
stimDiameterDeg = 16.6;   % degrees

% Generate stimulus template and set some defaults
if whichSite == 3
    experimentSpecs(3,1) = {'BAIR_ACER'};
end

stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);
stimParams.NIdaqRate = 1000;
stimParams.NIdaqNames = {'cDAQ1mod1'};

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
        
    case 'TACTILEVISUALSWEEP' %simple hand stimulation in different patterns for now
        
        directions = {'Descending'};
        conditions = {'TactileWeakVisualOffSweep', 'TactileStrongVisualOffSweep', 'TactileWeakVisualCongruentSweep', 'TactileWeakVisualIncongruentSweep'};
        stimParams.stimDurSecs            = 6; % seconds, length of stimulation at one finger
        stimParams.numSweeps              = 8; % how many sweeps across all fingers
        stimParams.interSweepIntervalSecs = 0; % pause in between sweeps in seconds
        stimParams.preScanDurSecs         = 12; % pause at the beginning of one run in secs
        stimParams.postScanDurSecs        = 12; % pause at the end of one run in secs
        stimParams.numOfStimulators       = 4;
        stimParams.modulatingFreq         = 6;
        stimParams.carrierFreq            = 30;
        
        %reset TR
        TR = 1;
        
        switch(lower(stimParams.modality))
            case 'fmri'
                stimParams.preScanDurSecs         = round(stimParams.preScanDurSecs/TR)*TR;
                stimParams.stimDurSecs            = round(stimParams.stimDurSecs/TR)*TR;
                stimParams.interSweepIntervalSecs = round(stimParams.interSweepIntervalSecs/TR)*TR;
            case {'ecog' 'eeg' 'meg'}
                stimParams.preScanDurSecs         = stimParams.preScanDurSecs; % seconds
            otherwise
                error('Unknown modality')
        end
        
        numberOfRuns = 4;
        
        %loop through the conditions
        for kk = 1:length(conditions)
            %set parameter accordingly
            stimParams.condition = conditions{kk};
            for runNum = 1:numberOfRuns
                % MAKE TASK EXPERIMENT
                stimMakeTactileVisualHandExperiment(stimParams, runNum,...
                    directions, makeFigure);
            end
        end
        
    case 'TACTILEVISUALBLOCKED' %simple hand stimulation in different patterns for now
        
        directions = {'Blocked'};
        conditions = {'TactileWeakVisualOffBlocked'};
        stimParams.stimDurSecs            = 12; % seconds, length of stimulation at one finger
        stimParams.numSweeps              = 8; % how many sweeps across all fingers
        stimParams.interSweepIntervalSecs = 12; % pause in between sweeps in seconds
        stimParams.preScanDurSecs         = 12; % pause at the beginning of one run in secs
        stimParams.postScanDurSecs        = 12;% pause at the end of one run in secs
        stimParams.numOfStimulators       = 4;
        stimParams.modulatingFreq         = 6;
        stimParams.carrierFreq            = 30;
        
        %reset TR
        TR = 1;
        
        switch(lower(stimParams.modality))
            case 'fmri'
                stimParams.preScanDurSecs         = round(stimParams.preScanDurSecs/TR)*TR;
                stimParams.stimDurSecs            = round(stimParams.stimDurSecs/TR)*TR;
                stimParams.interSweepIntervalSecs = round(stimParams.interSweepIntervalSecs/TR)*TR;
            case {'ecog' 'eeg' 'meg'}
                stimParams.preScanDurSecs         = stimParams.preScanDurSecs; % seconds
            otherwise
                error('Unknown modality')
        end
        
        numberOfRuns = 2;
        
        %loop through the conditions
        for kk = 1:length(conditions)
            %set parameter accordingly
            stimParams.condition = conditions{kk};
            for runNum = 1:numberOfRuns
                % MAKE TASK EXPERIMENT
                stimMakeTactileVisualHandExperiment(stimParams, runNum,...
                    directions, makeFigure);
            end
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



