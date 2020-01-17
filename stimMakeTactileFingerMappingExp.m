function stimMakeTactileFingerMappingExp(stimParams, runNum, directions, condition, makeFigure)
% Makes a simple tactile hand experiment where the stimulators vibrate in
% several patterns and assumes one stimulator per finger.
%

% set some stimulus properties related to the hardware
stimulus = [];
stimulus.NIdaqRate              = stimParams.NIdaqRate;
stimulus.NIdaqNames             = stimParams.NIdaqNames;
stimulus.numOfStimulators       = stimParams.numOfStimulators;
stimulus.frameRate              = stimParams.display.frameRate;
stimulus.cmap                   = stimParams.stimulus.cmap;
stimulus.srcRect                = stimParams.stimulus.srcRect;
stimulus.dstRect                = stimParams.stimulus.destRect;
stimulus.display                = stimParams.display;

stimParams.fingerIdx                = [1:stimulus.numOfStimulators, 6 7];
stimParams.fingers                  = {'thumb','index', 'middle', 'ring', 'little','all', 'blank'};
stimParams.stimDurSamples           = round(stimParams.stimDurSecs * stimParams.NIdaqRate);
stimParams.stimDurFrames            = stimParams.stimDurSecs * stimParams.display.frameRate;
stimParams.onDurSamples             = stimParams.onDurSecs * stimParams.NIdaqRate;
stimParams.offDurSamples            = stimParams.offDurSecs * stimParams.NIdaqRate;
stimParams.onOffReps                = stimParams.stimDurSamples  / (stimParams.onDurSamples + stimParams.offDurSamples);

stimParams.tactileIntensity         = 1.0; % maximal value is 2, check state of amplifier, too

switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % store how many samples needed for a trigger signal
        stimParams.triggerSamples = 50;
end

% check whether timings match
if (floor(stimParams.onOffReps) ~= stimParams.onOffReps)
    error('stimulus duration cannot be calculated');
end
% warn if samples don't match secs
if (stimParams.stimDurFrames / stimParams.display.frameRate ~= stimParams.stimDurSecs)
    error('stimulus duration cannot be calculated');
end

%% stimulus sequence timing
for jj = 1:length(directions)
    switch directions{jj}
        case 'Random'
            %duration of the whole experiment
            stimParams.expDurSecs = stimParams.numReps * ... %number of repetitions
                stimParams.numOfStimulators * ... %number of fingers
                (stimParams.stimDurSecs ...%duration of one stimulus
                + stimParams.interStimIntervalSecs)...%pause after each stimulus
                + stimParams.preScanDurSecs + stimParams.postScanDurSecs;%pauses at the beginning and the end of the experiment
            %durations of different (non-)stimulation periods
            stimParams.periodDursSecsVect  = [stimParams.preScanDurSecs, ... %pre-stimulation
                repmat(...
                [stimParams.stimDurSecs,... %one tactile stimulation
                stimParams.interStimIntervalSecs],... % pause after one stimulation
                1,stimParams.numReps * stimParams.numOfStimulators)... % repeat for each repetition and finger
                stimParams.postScanDurSecs];% add post-experimental pause
            %onsets of tactile stimuli(
            stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs, ... %first stimulus starts after pre-experiment pause
                repmat(...
                stimParams.stimDurSecs... %one stimulation
                + stimParams.interStimIntervalSecs,... % add pause after each sweep = stimulation
                1,(stimParams.numReps * stimParams.numOfStimulators)-1)... % repeat for all sweeps
                ]);
    end
end
% calculate in frames
stimParams.expDurFrames = round(stimParams.expDurSecs * stimParams.display.frameRate);%duration of the whole experiment in frames


stimParams.stimOnsetsSamples         = stimParams.stimOnsetsSecs*stimParams.NIdaqRate;%onsets of tactile stim in samples
stimParams.stimOnsetsSecsFrameNormed = round(stimParams.stimOnsetsSecs*stimParams.display.frameRate)/stimParams.display.frameRate;
stimParams.stimOnsetsFrames          = stimParams.stimOnsetsSecsFrameNormed*stimParams.display.frameRate;%onsets of new tactile stim in frames
stimParams.seqtiming                 = 0:(1/stimParams.display.frameRate):stimParams.expDurSecs;%frame onsets
stimParams.fixSeq                    = ones(size(stimParams.seqtiming)); %determines whether and which type of fixation is shown
stimParams.seq                       = ones(size(stimParams.seqtiming));
stimParams.trigSeq                   = zeros(size(stimParams.seqtiming));
stimParams.trigSeq(stimParams.stimOnsetsFrames + 1) = 1;
stimParams.trigSeq(1) = 1;
stimParams.trigSeq(end) = 1;

stimulus.seqtiming = stimParams.seqtiming;
stimulus.fixSeq = stimParams.fixSeq;
stimulus.seq = stimParams.seq;
stimulus.trigSeq = stimParams.trigSeq;

%debug (check timings of nidaq onsets vs psychtoolbox trigger onsets)
% figure;stem(stimulus.seqtiming, stimulus.trigSeq); hold on
% stem(stimParams.stimOnsetsSecs, ones(1,length(stimParams.stimOnsetsSecs))) 

% make a blank to insert between simulus presentations
screenRect          = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));
images              = ones([screenRect 1], 'uint8');
images(:,:,1)       = 127;
stimulus.images     = images;

%% Make stimulus for the experiment

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
%sample of the nidaq for one whole experiment
stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);


for jj = 1:length(directions)
    switch directions{jj}
        case 'Random'
            stimulatorOrder = repmat(1:stimParams.numOfStimulators,1, stimParams.numReps);
            stimulatorOrder = stimulatorOrder(randperm(length(stimulatorOrder)));
    end
end

for jj = 1:length(directions)
    
    % Loop through the stimulus sequence
    for ii = 1:length(stimParams.stimOnsetsSamples)
        % Insert the tactile signal at the respective time points and stimulator
        stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) + stimParams.stimDurSamples - 1, ...
            stimulatorOrder(ii)) = stimParams.stimulusSignal;
    end
    
    switch lower(stimParams.site)
        case 'nyuecog'
            % Write trigger sequence
            % add fake tactile channel
            stimParams.vibrotactileStimulus = [stimParams.vibrotactileStimulus, ...
                zeros(size(stimParams.vibrotactileStimulus,1),1)];
            % Insert the trigger signal at the onset of the run
            stimParams.vibrotactileStimulus(int64(1):int64(1 +stimParams.triggerSamples - 1), ...
                stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
            % Loop through the stimulus sequence
            for ii = 1:length(stimParams.stimOnsetsSamples)
                % Insert the trigger signal at the respective time points
                % of the stimulator
                stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) +stimParams.triggerSamples - 1, ...
                    stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
            end
            % Insert the trigger signal at the end of the run
            stimParams.vibrotactileStimulus(end - stimParams.triggerSamples : end - 1, ...
                stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
    end
    
    %Save stimulus matrix, tsv info and make figures
    fname = sprintf('%s_%s%s_%d.mat', stimParams.site, condition, directions{jj}, runNum);
    
    % check whether figures should be made
    if makeFigure
        % save figure with images of stimulus
        f = figure('visible', 'off');
        imagesc(stimParams.vibrotactileStimulus')
        title (sprintf('%s', directions{jj}))
        ylabel('Tactile stimulators');
        xlabel('Time (s)')
        yticks(stimParams.fingerIdx(1:stimParams.numOfStimulators))
        yticklabels(stimParams.fingers(1:stimParams.numOfStimulators))
        xticks(0:stimParams.NIdaqRate * 10:stimParams.expDurSecs * stimParams.NIdaqRate)
        xticklabels(0:10:stimParams.expDurSecs)
        
        saveas(f, fullfile(vistadispRootPath, 'StimFiles', sprintf('%s.png',fname(1:end-6))))
    end
    
    % Sparsify the visual stimulus sequence and the triggers (apart
    % from the tactile one)
    maxUpdateInterval = 0.25;
    stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);
    
    % Set TSV file information
    onsets          = round(stimParams.stimOnsetsSecs,3)';
    duration        = repmat(stimParams.stimDurSecs, length(stimParams.stimOnsetsSecs),1);
    ISI             = repmat(stimParams.offDurSecs, length(stimParams.stimOnsetsSecs),1);
    trial_type      = stimParams.fingerIdx(stimulatorOrder)';
    trial_name      = stimParams.fingers(stimulatorOrder)';
    stim_frequency  = repmat(stimParams.carrierFreq, length(stimParams.stimOnsetsSecs),1);
    stim_amplitude  = ones(length(stimParams.stimOnsetsSecs),1);
    stim_file       = repmat(fname, length(stimParams.stimOnsetsSecs),1);
    stim_order      = repmat(directions(jj), length(stimParams.stimOnsetsSecs),1);
    
    stimParams.tsv = table(onsets, duration, ISI, trial_type, trial_name, ...
        stim_frequency, stim_amplitude, stim_file, stim_order);
    stimulus.tsv = table(onsets, duration, ISI, trial_type, trial_name, ...
        stim_frequency, stim_amplitude, stim_file, stim_order);
    
        
    % store for main script calls
    stimulus.vibrotactileStimulus = stimParams.vibrotactileStimulus;
    
    % Save
    if ~exist(fullfile(vistadispRootPath, 'StimFiles',  fname), 'file')
        fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
        save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')
        save(fullfile(vistadispRootPath, 'StimFiles',  [sprintf('%s_%s', stimParams.site, directions{jj}(1:4)),'_stimParams.mat']), 'stimParams');
    else
        error('stimulus files for this experiment already exist, move or delete old files')
    end
end

