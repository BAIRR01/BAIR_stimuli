function stimMakeTactileTemporalExp(stimParams, runNum, condition, makeFigure)
% Stimulators vibrate for a pre-defined duration or two vibrations are separated
% by the pre-defined duration


switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % store how many samples are needed for a trigger signal
        stimParams.triggerSamples = 50;
end

% calculate number of samples for each stimulus (in terms of the NIdaq device)
stimParams.testedDurSamples         = round(stimParams.testedDurSecs * stimParams.NIdaqRate);
stimParams.tapDurSamples            = stimParams.tapDurSecs * stimParams.NIdaqRate;


% create condition matrix
stimParams.conditionsUnique = combvec(stimParams.testedDurSecs, stimParams.tapCondition)';
% add condition identifier
stimParams.conditionsUnique = [stimParams.conditionsUnique, (1:size(stimParams.conditionsUnique,1))'];


%% stimulation sequence
% read stimulus conditions order and onsets from file
stimParams.stimSeqOptmized = dlmread(fullfile(vistadispRootPath, 'StimFiles',  sprintf('tactileTemporalSeq_%d.txt', runNum)));


switch(lower(stimParams.modality))
    case 'fmri'
        % read out condition info
        stimParams.conditions = stimParams.stimSeqOptmized(:,[3,4,2]);
        % read out stimulus onsets
        stimParams.stimOnsetsSecs = stimParams.preScanDurSecs + stimParams.stimSeqOptmized(:,1);
        % adjust for TR
        stimParams.stimOnsetsSecs = round(stimParams.stimOnsetsSecs/(0.2*stimParams.TR))*(0.2*stimParams.TR);
        % calculate ITI (subtract stimulus duration + 2 * tap duration in
        % two tap condition)
        stimParams.ITIs = diff(stimParams.stimOnsetsSecs)-stimParams.conditions(1:end-1,1)-(stimParams.conditions(1:end-1,2)-1)*2*stimParams.tapDurSecs;
    case {'ecog' 'eeg' 'meg'}
        % random order of stimuli
        stimParams.conditions = repmat(stimParams.conditionsUnique, stimParams.numReps,1);
        stimParams.conditions = stimParams.conditions(randperm(size(stimParams.conditions, 1)), :);
        % random ITIs
        stimParams.ITIs = linspace(stimParams.minITI,stimParams.maxITI,size(stimParams.conditions, 1)-1);
        stimParams.ITIs = stimParams.ITIs(randperm(length(stimParams.ITIs)));
        %onsets of tactile stimuli(
        stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs; ... % first stimulus starts after pre-experiment pause
            stimParams.conditions(1:end-1,1) + (stimParams.conditions(1:end-1,2) == 2) * 2 * stimParams.tapDurSecs ... % stimulus durations
            + stimParams.ITIs'... % ITI
            ]);
    otherwise
        error('Unknown modality')
end

%duration of the whole experiment
stimParams.expDurSecs = sum(stimParams.conditions(:, 1)) ... % durations of the stimulus(either continuous vibration or gap between 2 taps)
    + length(stimParams.conditions(stimParams.conditions(:,2) == 2)) * 2 * stimParams.tapDurSecs ... % duration of the taps in 2 tap conditions
    + sum(stimParams.ITIs) ...% pause after each stimulus
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

%initialize matrix to hold activation of each stimulator during each
%sample of the nidaq for one whole experiment
stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);

% Loop through the stimulus sequence
for ii = 1:length(stimParams.stimOnsetsSamples)
    % check whether this is 1 tap or 2 taps stimulus
    if stimParams.conditions(ii,2) == 1
        % Insert the tactile signal at the respective time points
        stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):...
            int64(stimParams.stimOnsetsSamples(ii) + 1) + stimParams.conditions(ii,1) * stimParams.NIdaqRate - 1, ...
            :) = ...
            repmat(0.5 + 0.5 * sin(-pi/2+...
            linspace(0, ...
            (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.conditions(ii,1) * stimParams.NIdaqRate), ...
            stimParams.conditions(ii,1) * stimParams.NIdaqRate)'...
            ), ...
            1,stimParams.numOfStimulators);
    elseif stimParams.conditions(ii,2) == 2
        % Insert the tactile signal for the first tap
        stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):...
            int64(stimParams.stimOnsetsSamples(ii) + 1 + stimParams.tapDurSamples - 1), ...
            :) = repmat(0.5 + 0.5 * sin(-pi/2+linspace(0, ...
            (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.tapDurSamples), ...
            stimParams.tapDurSamples)'),1,stimParams.numOfStimulators);
        % Insert the tactile signal for the second tap
        stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1 ...
            + stimParams.tapDurSamples + stimParams.conditions(ii,1) * stimParams.NIdaqRate):...
            int64(stimParams.stimOnsetsSamples(ii) + 1 + stimParams.tapDurSamples ...
            + stimParams.conditions(ii,1) * stimParams.NIdaqRate ...
            + stimParams.tapDurSamples - 1), ...
            :) = repmat(0.5 + 0.5 * sin(-pi/2+linspace(0, ...
            (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.tapDurSamples), ...
            stimParams.tapDurSamples)'),1,stimParams.numOfStimulators);
    end
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
duration        = stimParams.conditions(:,1) .* (stimParams.conditions(:,2) == 1) ...
    + stimParams.tapDurSecs * (stimParams.conditions(:,2) == 2);
ISI             = stimParams.conditions(:,1) .* (stimParams.conditions(:,2) == 2);
trial_type      = stimParams.conditionNumbers(stimParams.conditions(:,3))';
trial_name      = stimParams.conditionNames(stimParams.conditions(:,3))';
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

