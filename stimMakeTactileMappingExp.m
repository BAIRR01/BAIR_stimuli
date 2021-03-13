function stimMakeTactileMappingExp(stimParams, runNum, condition, makeFigure)
% Stimulators in pre-specified locations vibrate in random order


switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % store how many samples are needed for a trigger signal
        stimParams.triggerSamples = 50;
end

% calculate number of samples for each stimulus (in terms of the NIdaq device)
stimParams.stimDurSecs      = (stimParams.onDurSecs + stimParams.offDurSecs) * stimParams.onOffReps;
stimParams.stimDurSamples   = round(stimParams.stimDurSecs * stimParams.NIdaqRate);
stimParams.stimDurFrames    = stimParams.stimDurSecs * stimParams.display.frameRate;
stimParams.onDurSamples     = stimParams.onDurSecs * stimParams.NIdaqRate;
stimParams.offDurSamples    = stimParams.offDurSecs * stimParams.NIdaqRate;


%% stimulation sequence
% read stimulus conditions order and onsets from file
stimParams.stimSeqOptmized = dlmread(fullfile(vistadispRootPath, 'StimFiles',  sprintf('tactileSpatialSeq_%d.txt', runNum)));
% fixed ITIs
stimParams.ITIs = repmat(stimParams.ITI,length(stimParams.fingerIdx)* stimParams.numReps-1);

% determine stimulus onsets (depends on modality)
switch(lower(stimParams.modality))
    case 'fmri'
        % read out condition info
        stimParams.conditions = stimParams.stimSeqOptmized(:,2);
        % read out stimulus onsets
        stimParams.stimOnsetsSecs = stimParams.preScanDurSecs + stimParams.stimSeqOptmized(:,1);
        % adjust for TR
        stimParams.stimOnsetsSecs = round(stimParams.stimOnsetsSecs/(0.2*stimParams.TR))*(0.2*stimParams.TR);
    case {'ecog' 'eeg' 'meg'}
        % random order stimuli
        stimParams.conditions = repmat(stimParams.fingerIdx,1, stimParams.numReps);
        stimParams.conditions = stimParams.conditions(randperm(size(stimParams.conditions, 2)))';
        % onsets of tactile stimuli(
        stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs; ... % first stimulus starts after pre-experiment pause
            repmat(...
            stimParams.stimDurSecs... %one stimulation
            + stimParams.ITI,... % add ITI
            size(stimParams.conditions,1)-1,1).... % repeat for all stimuli
            ]);
    otherwise
        error('Unknown modality')
end

% duration of the whole experiment
stimParams.expDurSecs = stimParams.numReps * ... % number of repetitions
    stimParams.numOfStimulators * ... % number of stimulators
    (stimParams.stimDurSecs ...% duration of one stimulus
    + stimParams.ITI)...% pause after each stimulus
    + stimParams.preScanDurSecs + stimParams.postScanDurSecs;% pauses at the beginning and the end of the experiment


% calculate experiment duration in frames
stimParams.expDurFrames = round(stimParams.expDurSecs * stimParams.display.frameRate);%duration of the whole experiment in frames

% convert units
stimParams.stimOnsetsSamples         = stimParams.stimOnsetsSecs*stimParams.NIdaqRate;%onsets of tactile stim in samples
stimParams.stimOnsetsSecsFrameNormed = round(stimParams.stimOnsetsSecs*stimParams.display.frameRate)/stimParams.display.frameRate;
stimParams.stimOnsetsFrames          = round(stimParams.stimOnsetsSecsFrameNormed*stimParams.display.frameRate);%onsets of new tactile stim in frames
% create timing variable analogous to the visual experiments (main loop
% runs on frames)
stimParams.seqtiming                 = 0:(1/stimParams.display.frameRate):stimParams.expDurSecs;%frame onsets
stimParams.fixSeq                    = ones(size(stimParams.seqtiming)); %determines whether and which type of fixation is shown
stimParams.seq                       = ones(size(stimParams.seqtiming));
% add trigger timing (visual, auditory)
stimParams.trigSeq                   = zeros(size(stimParams.seqtiming));
stimParams.trigSeq(stimParams.stimOnsetsFrames + 1) = 1;
stimParams.trigSeq(1) = 1;
stimParams.trigSeq(end) = 1;

% make a blank to insert between simulus presentations
screenRect          = size(zeros(stimParams.stimulus.destRect(4)-stimParams.stimulus.destRect(2): stimParams.stimulus.destRect(3)-stimParams.stimulus.destRect(1)));
images              = ones([screenRect 1], 'uint8');
images(:,:,1)       = 127;
stimParams.images   = images;


%% make stimulus
% signal for one tactile stimulus (e.g., per finger)
% carrier signal for one tactile stimulus
stimulusSignalBase = 0.5 + 0.5 * sin(-pi/2+linspace(0, ...
    (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.stimDurSamples), ...
    stimParams.stimDurSamples)');
% modulating signal for one tactile stimulus
stimulusSignalModulator = repmat([ones(stimParams.onDurSamples, 1); zeros(stimParams.offDurSamples,1)], ...
    stimParams.onOffReps,1);
% final signal for one tactile stimulus
stimParams.stimulusSignal = stimulusSignalBase .* stimulusSignalModulator * stimParams.tactileIntensity;
% replace last entry with 0 to inactivate tactile stimulator
stimParams.stimulusSignal(end) = 0;

%initialize matrix to hold activation of each stimulator during each
%sample of the nidaq for one run
stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);

%initialize matrix to hold activation of each stimulator during each
%sample of the nidaq for one whole experiment
stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);

% Loop through the stimulus sequence
for ii = 1:length(stimParams.stimOnsetsSamples)
    % Insert the tactile signal at the respective time points and stimulator
    stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) + stimParams.stimDurSamples - 1, ...
        stimParams.conditions(ii)) = stimParams.stimulusSignal;
end

switch lower(stimParams.site)
    case 'nyuecog'
        % Write trigger sequence
        % add fake tactile signal
        stimParams.vibrotactileStimulus = [stimParams.vibrotactileStimulus, ...
            zeros(size(stimParams.vibrotactileStimulus,1),1)];
        % Insert the trigger signal at the onset of the run
        stimParams.vibrotactileStimulus(int64(1):int64(1 +stimParams.triggerSamples - 1), ...
            stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
        % Loop through the stimulus sequence
        for ii = 1:length(stimParams.stimOnsetsSamples)
            % Insert the trigger signal at the onset of each stimulus (not
            % of each tap)
            stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) +stimParams.triggerSamples - 1, ...
                stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
        end
        % Insert the trigger signal at the end of the run
        stimParams.vibrotactileStimulus(end - stimParams.triggerSamples : end - 1, ...
            stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
end

%% save stimulus matrix, tsv info, and make figures

% name of the stimulus file
fname = sprintf('%s_%s_%d.mat', stimParams.site, condition, runNum);

% check whether figures should be made
if makeFigure
    % save figure with images of stimulus
    f = figure('visible', 'off');
    imagesc(stimParams.vibrotactileStimulus'>0)
    title (sprintf('%s', condition))
    ylabel('Tactile stimulators');
    xlabel('Time (s)')
    yticks(stimParams.fingerIdx(1:stimParams.numOfStimulators))
    yticklabels(stimParams.fingers(1:stimParams.numOfStimulators))
    xticks(0:stimParams.NIdaqRate * 10:stimParams.expDurSecs * stimParams.NIdaqRate)
    xticklabels(0:10:stimParams.expDurSecs)
    set(gcf,'units','centimeters','position',[0,0,50,10])
    
    saveas(f, fullfile(vistadispRootPath, 'StimFiles', sprintf('%s.png',fname(1:end-6))))
end

% Sparsify the visual stimulus sequence and the triggers (apart
% from the tactile one)
maxUpdateInterval = 0.25;
stimParams = sparsifyStimulusStruct(stimParams, maxUpdateInterval);

% Set TSV file information
onsets          = round(stimParams.stimOnsetsSecs,3);
duration        = repmat(stimParams.stimDurSecs, size(stimParams.stimOnsetsSecs,1), 1);
ISI             = repmat(stimParams.ITI, size(stimParams.stimOnsetsSecs,1), 1);
trial_type      = stimParams.conditionNumbers(stimParams.conditions)';
trial_name      = stimParams.conditionNames(stimParams.conditions)';
stim_frequency  = repmat(stimParams.carrierFreq, length(stimParams.stimOnsetsSecs),1);
stim_amplitude  = ones(length(stimParams.stimOnsetsSecs),1);
stim_file       = repmat({fname}, length(stimParams.stimOnsetsSecs),1);

stimParams.tsv = table(onsets, duration, ISI, trial_type, trial_name, ...
    stim_frequency, stim_amplitude, stim_file);

% store all information needed to run the experiment (i.e., information the
% experimental script is loading) in a structure called 'stimulus'
stimulus = [];
stimulus.NIdaqRate              = stimParams.NIdaqRate;
stimulus.NIdaqNames             = stimParams.NIdaqNames;
stimulus.numOfStimulators       = stimParams.numOfStimulators;
stimulus.frameRate              = stimParams.display.frameRate;
stimulus.cmap                   = stimParams.stimulus.cmap;
stimulus.srcRect                = stimParams.stimulus.srcRect;
stimulus.dstRect                = stimParams.stimulus.destRect;
stimulus.display                = stimParams.display;
stimulus.seqtiming              = stimParams.seqtiming;
stimulus.fixSeq                 = stimParams.fixSeq;
stimulus.seq                    = stimParams.seq;
stimulus.trigSeq                = stimParams.trigSeq;
stimulus.images                 = stimParams.images;

stimulus.tsv = table(onsets, duration, ISI, trial_type, trial_name, ...
    stim_frequency, stim_amplitude, stim_file);

stimulus.vibrotactileStimulus = stimParams.vibrotactileStimulus;

% save the stimulus structure to be loaded whenever the experiment is run
% (note: this file is not synchronized via git, it has to be created on
% each computer that runs the experiment)
if ~exist(fullfile(vistadispRootPath, 'StimFiles',  fname), 'file')
    fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
    save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')
else
    error('stimulus files for this experiment already exist, move or delete old files')
end

