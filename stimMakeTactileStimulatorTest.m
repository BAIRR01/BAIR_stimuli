function stimMakeTactileStimulatorTest(stimParams, runNum, directions, condition, makeFigure)
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

stimParams.fingerIdx                = [1:stimulus.numOfStimulators];
stimParams.fingers                  = {'thumb','index', 'middle', 'ring', 'little'};

stimParams.tactileIntensity         = 1.0; % maximal value is 2, check state of amplifier, too


switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % store how many samples are needed for a trigger signal
        stimParams.triggerSamples = 50;
end

stimParams.stimDurSamples            = stimParams.stimDurSecs * stimParams.NIdaqRate;


%% stimulus sequence timing
for jj = 1:length(directions)
    switch directions{jj}
        case {'Ascending', 'Descending'}
            %% stimulus sequence timing
            %duration of the whole experiment
            stimParams.expDurSecs = stimParams.numReps * ... %number of repetitions
                stimParams.numOfStimulators * ... %number of fingers
                (stimParams.stimDurSecs ...%duration of one stimulus
                + stimParams.interStimIntervalSecs)...%pause after each stimulus
                + stimParams.preScanDurSecs + stimParams.postScanDurSecs;%pauses at the beginning and the end of the experiment
            
            %onsets of tactile stimuli(
            stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs, ... %first stimulus starts after pre-experiment pause
                repmat(...
                stimParams.stimDurSecs... %one stimulation
                + stimParams.interStimIntervalSecs,... % add pause after each sweep = stimulation
                1,(stimParams.numReps * stimParams.numOfStimulators)-1)... % repeat for all sweeps
                ]);
        case 'All'
            %% stimulus sequence timing
            %duration of the whole experiment
            stimParams.expDurSecs = stimParams.numReps * ... %number of repetitions
                (stimParams.stimDurSecs ...%duration of one stimulus
                + stimParams.interStimIntervalSecs)...%pause after each stimulus
                + stimParams.preScanDurSecs + stimParams.postScanDurSecs;%pauses at the beginning and the end of the experiment
            
            %onsets of tactile stimuli(
            stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs, ... %first stimulus starts after pre-experiment pause
                repmat(...
                stimParams.stimDurSecs... %one stimulation
                + stimParams.interStimIntervalSecs,... % add pause after each sweep = stimulation
                1,(stimParams.numReps)-1)... % repeat for all sweeps
                ]);
    end
    
    % calculate experiment duration in frames
    stimParams.expDurFrames = round(stimParams.expDurSecs * stimParams.display.frameRate);%duration of the whole experiment in frames
    
    % convert units
    stimParams.stimOnsetsSamples         = stimParams.stimOnsetsSecs*stimParams.NIdaqRate;%onsets of tactile stim in samples
    stimParams.stimOnsetsSecsFrameNormed = round(stimParams.stimOnsetsSecs*stimParams.display.frameRate)/stimParams.display.frameRate;
    stimParams.stimOnsetsFrames          = round(stimParams.stimOnsetsSecsFrameNormed*stimParams.display.frameRate);%onsets of new tactile stim in frames
    % create timing analogous to the visual one
    stimParams.seqtiming                 = 0:(1/stimParams.display.frameRate):stimParams.expDurSecs;%frame onsets
    stimParams.fixSeq                    = ones(size(stimParams.seqtiming)); %determines whether and which type of fixation is shown
    stimParams.seq                       = ones(size(stimParams.seqtiming));
    % add trigger timing (visual, auditory)
    stimParams.trigSeq                   = zeros(size(stimParams.seqtiming));
    stimParams.trigSeq(stimParams.stimOnsetsFrames + 1) = 1;
    stimParams.trigSeq(1) = 1;
    stimParams.trigSeq(end) = 1;
    % store in stimulus structure
    stimulus.seqtiming = stimParams.seqtiming;
    stimulus.fixSeq = stimParams.fixSeq;
    stimulus.seq = stimParams.seq;
    stimulus.trigSeq = stimParams.trigSeq;
    
    % make a blank to insert between simulus presentations
    screenRect          = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));
    images              = ones([screenRect 1], 'uint8');
    images(:,:,1)       = 127;
    stimulus.images     = images;
    
    %% Make stimulus for the experiment
    
    %initialize matrix to hold activation of each stimulator during each
    %sample of the nidaq for one whole experiment
    stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);
    
    % signal for one tactile stimulus (e.g., per finger)
    % carrier signal for one tactile stimulus
    stimParams.stimulusSignal = stimParams.tactileIntensity * (0.5 + 0.5 * sin(-pi/2+linspace(0, ...
        (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.stimDurSamples), ...
        stimParams.stimDurSamples)'));
    % replace last entry with 0 to inactivate tactile stimulator
    stimParams.stimulusSignal(end) = 0;
    
    switch directions{jj}
        case 'Ascending'
            stimulatorOrder = repmat(1:stimParams.numOfStimulators,1, stimParams.numReps);
        case 'Descending'
            stimulatorOrder = repmat(stimParams.numOfStimulators:-1:1,1, stimParams.numReps);
    end
    
    switch directions{jj}
        case {'Ascending', 'Descending'}
            % Loop through the stimulus sequence
            for ii = 1:length(stimParams.stimOnsetsSamples)
                % Insert the tactile signal at the respective time points and stimulator
                stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) + stimParams.stimDurSamples - 1, ...
                    stimulatorOrder(ii)) = stimParams.stimulusSignal;
            end
        case 'All'
            for ii = 1:length(stimParams.stimOnsetsSamples)
                % Insert the tactile signal at the respective time points and stimulator
                stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) + stimParams.stimDurSamples - 1, ...
                    :) = repmat(stimParams.stimulusSignal, 1, stimParams.numOfStimulators);
            end
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
                % Insert the tactile signal at the respective time points and stimulator
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
        title (sprintf('%s', condition))
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
    ISI             = repmat(stimParams.interStimIntervalSecs, length(stimParams.stimOnsetsSecs),1);
    switch directions{jj}
        case {'Ascending', 'Descending'}
            trial_type      = stimParams.fingerIdx(stimulatorOrder)';
            trial_name      = stimParams.fingers(stimulatorOrder)';
        case 'All'
            trial_type      = repmat(100,length(stimParams.stimOnsetsSecs),1);
            trial_name      = repmat('All',length(stimParams.stimOnsetsSecs),1);
    end
    stim_frequency  = repmat(stimParams.carrierFreq, length(stimParams.stimOnsetsSecs),1);
    stim_amplitude  = repmat(stimParams.tactileIntensity,length(stimParams.stimOnsetsSecs),1);
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
    else
        error('stimulus files for this experiment already exist, move or delete old files')
    end
    
end
