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
tmp = daq.getDevices;
stimParams.NIdaqNames = tmp.ID;

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
        conditions = {'TactileWeakVisualOffSweep', 'TactileStrongVisualOffSweep',...
            'TactileWeakVisualCongruentSweep', 'TactileWeakVisualShiftToLittleFingerSweep', 'TactileWeakVisualShiftToThumbSweep'};
        stimParams.stimDurSecs            = 4; % seconds, length of stimulation at one finger
        stimParams.numSweeps              = 6; % how many sweeps across all fingers
        stimParams.interSweepIntervalSecs = 0; % pause in between sweeps in seconds
        stimParams.preScanDurSecs         = 0; % pause at the beginning of one run in secs
        stimParams.postScanDurSecs        = 10; % pause at the end of one run in secs
        stimParams.numOfStimulators       = 5;
        stimParams.modulatingFreq         = 4;
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
        
        numberOfRuns = 8;
        
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
        
        
    case 'RANDOM_5FINGERS' % single fingers in random order
        
        directions = {'Random'};
        condition = 'TactileFingerMapping';
        stimParams.stimDurSecs              = 1; % seconds, length of stimulation at one finger
        stimParams.numReps                  = 10; % how many repetitions across all fingers in one run
        stimParams.interStimIntervalSecs    = 3; % pause in between single stimuli in seconds
        stimParams.preScanDurSecs           = 3; % pause at the beginning of one run in secs
        stimParams.postScanDurSecs          = 3; % pause at the end of one run in secs
        stimParams.numOfStimulators         = 5;
        stimParams.carrierFreq              = 110; % base vibration in Hz
        stimParams.onDurSecs                = 0.4; % duration of constant vibration
        stimParams.offDurSecs               = 0.1; % duration of break betweeen constant vibrations
        
        
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
        
        for runNum = 1:numberOfRuns
            % make stimulus for experiment
            stimMakeTactileFingerMappingExp(stimParams, runNum,...
                directions, condition, makeFigure);
        end
        
    case 'TEMPORAL' % all fingers, different durations
        
        condition = 'TactileTemporal';
        stimParams.testedDurSecs            = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2]; % seconds, length of tested stimulation (either constant vibration or gap between vibrations)
        stimParams.tapCondition             = [1, 2]; % either constant vibration == 1 or gap between two vibrations == 2
        stimParams.tapDurSecs               = 0.2; % seconds, duration of the taps in 2 tap condition
        
        stimParams.conditionNames = {...
            'ONE-PULSE-1','ONE-PULSE-2','ONE-PULSE-3', ...
            'ONE-PULSE-4','ONE-PULSE-5','ONE-PULSE-6',...
            'TWO-PULSE-1','TWO-PULSE-2','TWO-PULSE-3', ...
            'TWO-PULSE-4','TWO-PULSE-5','TWO-PULSE-6'};
        stimParams.conditionNumbers = 119:130;
        
        stimParams.carrierFreq              = 110; % base vibration in Hz
        stimParams.tactileIntensity         = 1.0; % maximal value is 2, check state of amplifier, too
        
        
        switch(lower(stimParams.modality))
            case 'fmri'
                stimParams.numReps              = 6; % how many repetitions across all fingers in one run
                numberOfRuns                    = 2;
                stimParams.preScanDurSecs       = 11.9;
                stimParams.postScanDurSecs      = 11.9;
                stimParams.TR                   = 0.85;
            case {'ecog' 'eeg' 'meg'}
                stimParams.numReps              = 12; % how many repetitions across all fingers in one run
                numberOfRuns                    = 1;
                stimParams.preScanDurSecs       = 3; % pause at the beginning of one run in secs
                stimParams.postScanDurSecs      = 3; % pause at the end of one run in secs
                stimParams.minITI               = 1.25; % interval between stimulus offset and next stimulus onset
                stimParams.maxITI               = 1.75; % interval between stimulus offset and next stimulus onset
            otherwise
                error('Unknown modality')
        end
        
        
        stimParams.numOfStimulators         = 5; % no of stimulators (at NYU ECoG one additional tactile channel serves as trigger channel)
        stimParams.fingerIdx                = 1:stimParams.numOfStimulators;
        stimParams.fingers                  = {'thumb','index', 'middle', 'ring', 'little'}; % name the stimulators (make sure the names, locations, and stimulators match)
        
        for runNum = 1:numberOfRuns
            % make stimulus for experiment
            stimMakeTactileTemporalExp(stimParams, runNum,...
                condition, makeFigure);
        end
        
        
    case 'TACTILEMAPPING' % different fingers (and other locations) in a random design
        
        condition = 'TactileFingerMapping';
        stimParams.numOfStimulators         = 8; % no of stimulators (at NYU ECoG one additional tactile channel serves as trigger channel)
        stimParams.fingerIdx                = 1:stimParams.numOfStimulators;
        stimParams.fingers                  = {...
            'thumb','index', 'middle',...
            'ring', 'little', ...
            'index_dorsal_2', 'index_dorsal_3', 'little_dorsal_2'}; % name the stimulators (make sure the names, locations, and stimulators match)
        
        stimParams.conditionNames           = {...
            'thumb','index', 'middle',...
            'ring', 'little', ...
            'index_dorsal_2', 'index_dorsal_3', 'little_dorsal_2'};
        stimParams.conditionNumbers = 130+(1:length(stimParams.conditionNames));
        
        
        stimParams.carrierFreq              = 110; % base vibration in Hz
        stimParams.tactileIntensity         = 1.0; % maximal value is 2, check state of amplifier, too
        
        
        switch(lower(stimParams.modality))
            case 'fmri'
                stimParams.numReps              = 10; % how many repetitions of all locations in one run
                numberOfRuns                    = 2;
                stimParams.preScanDurSecs       = 2.55;
                stimParams.postScanDurSecs      = 2.55;
                stimParams.ITI                  = 2.55; % interval between stimulus offset and next stimulus onset
                stimParams.onDurSecs            = 0.350; % duration of constant vibration
                stimParams.offDurSecs           = 0.075; % duration of break betweeen constant vibrations
                stimParams.onOffReps            = 6; % no repetitions of onDur + off Dur
                stimParams.TR                   = 0.85;
            case {'ecog' 'eeg' 'meg'}
                stimParams.numReps              = 20; % how many repetitions of all locations in one run
                numberOfRuns                    = 1;
                stimParams.preScanDurSecs       = 3; % pause at the beginning of one run in secs
                stimParams.postScanDurSecs      = 3; % pause at the end of one run in secs
                stimParams.ITI                  = 1; % interval between stimulus offset and next stimulus onset
                stimParams.onDurSecs            = 0.350; % duration of constant vibration
                stimParams.offDurSecs           = 0.75; % duration of break betweeen constant vibrations
                stimParams.onOffReps            = 2; % no repetitions of onDur + off Dur
            otherwise
                error('Unknown modality')
        end
        
        
        for runNum = 1:numberOfRuns
            % make stimulus for experiment
            stimMakeTactileMappingExp(stimParams, runNum,...
                condition, makeFigure);
        end
        
        
        
    case 'STIMULATORTEST' % all fingers, different durations
        
        directions = {'All', 'Ascending', 'Descending'};
        condition = 'TactTest';
        
        stimParams.stimDurSecs              = 1; % seconds, length of stimulation at one location
        stimParams.interStimIntervalSecs    = 0.5; % pause in between single stimuli in seconds
        stimParams.numReps                  = 2; % how many repetitions across all fingers in one run
        stimParams.preScanDurSecs           = 2; % pause at the beginning of one run in secs
        stimParams.postScanDurSecs          = 2; % pause at the end of one run in secs
        stimParams.numOfStimulators         = 8;
        stimParams.fingerIdx                = 1:stimParams.numOfStimulators;
        stimParams.fingers                  = {...
            'thumb','index', 'middle',...
            'ring', 'little', ...
            'index_dorsal_2', 'index_dorsal_3', 'little_dorsal_2'}; % name the stimulators (make sure the names, locations, and stimulators match)
        
        stimParams.conditionNames           = {...
            'thumb','index', 'middle',...
            'ring', 'little', ...
            'index_dorsal_2', 'index_dorsal_3', 'little_dorsal_2'};
        
        stimParams.carrierFreq              = 110; % base vibration in Hz
        stimParams.tactileIntensity         = 1.0; % maximal value is 2, check state of amplifier, too
        
        
        numberOfRuns = 1;
        
        for runNum = 1:numberOfRuns
            % make stimulus for experiment
            stimMakeTactileStimulatorTest(stimParams, runNum,...
                directions, condition, makeFigure);
        end
        
end



