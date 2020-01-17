function stimMakeTactileTemporalExp(stimParams, runNum, condition, makeFigure)
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


stimParams.testedDurSamples         = round(stimParams.testedDurSecs * stimParams.NIdaqRate);
stimParams.tapDurSamples            = stimParams.tapDurSecs * stimParams.NIdaqRate;

% create condition matrix
stimParams.conditions = repmat(combvec(stimParams.testedDurSecs, stimParams.tapCondition),1, stimParams.numReps)';
% shuffle conditions
stimParams.conditions = stimParams.conditions(randperm(size(stimParams.conditions, 1)), :);

% random ISIs
stimParams.ISIs = datasample(1+[1 2 4 8 16 32]/stimParams.display.frameRate, ...
    size(stimParams.conditions, 1)-1);


%% stimulus sequence timing
%duration of the whole experiment
stimParams.expDurSecs = sum(stimParams.conditions(:, 1)) ... %durations (1 and 2 tap)
    + length(stimParams.conditions(stimParams.conditions(:,2) == 2)) * 2 * stimParams.tapDurSecs ... % duration of the taps in 2 tap conditions
    + sum(stimParams.ISIs) ...%pause after each stimulus
    + stimParams.preScanDurSecs + stimParams.postScanDurSecs;%pauses at the beginning and the end of the experiment
%onsets of tactile stimuli(
stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs; ... %first stimulus starts after pre-experiment pause
    stimParams.conditions(1:end-1,1) + (stimParams.conditions(1:end-1,2) == 2) * 2 * stimParams.tapDurSecs ... % stimulus durations
    + [0,stimParams.ISIs(1:end-1)]'... % ISI
    ]);

% calculate in frames
stimParams.expDurFrames = round(stimParams.expDurSecs * stimParams.display.frameRate);%duration of the whole experiment in frames


stimParams.stimOnsetsSamples         = stimParams.stimOnsetsSecs*stimParams.NIdaqRate;%onsets of tactile stim in samples
stimParams.stimOnsetsSecsFrameNormed = round(stimParams.stimOnsetsSecs*stimParams.display.frameRate)/stimParams.display.frameRate;
stimParams.stimOnsetsFrames          = stimParams.stimOnsetsSecsFrameNormed*stimParams.display.frameRate;%onsets of new tactile stim in frames
stimParams.seqtiming                 = 0:(1/stimParams.display.frameRate):stimParams.expDurSecs;%frame onsets
stimParams.fixSeq                    = ones(size(stimParams.seqtiming)); %determines whether and which type of fixation is shown
stimParams.seq                       = ones(size(stimParams.seqtiming));

stimulus.seqtiming = stimParams.seqtiming;
stimulus.fixSeq = stimParams.fixSeq;
stimulus.seq = stimParams.seq;

% make a blank to insert between simulus presentations
screenRect          = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));
images              = ones([screenRect 1], 'uint8');
images(:,:,1)       = 127;
stimulus.images     = images;

%% Make stimulus for the experiment

%initialize matrix to hold activation of each stimulator during each
%sample of the nidaq for one whole experiment
stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);

% Loop through the stimulus sequence
for ii = 1:length(stimParams.stimOnsetsSamples)
    % check whether this is 1 tap or 2 taps
    if stimParams.conditions(ii,2) == 1
        % Insert the tactile signal at the respective time points
        stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):...
            int64(stimParams.stimOnsetsSamples(ii) + 1) + stimParams.conditions(ii,1) * stimParams.NIdaqRate - 1, ...
            :) = repmat(0.5 + 0.5 * sin(-pi/2+linspace(0, ...
            (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.conditions(ii,1) * stimParams.NIdaqRate), ...
            stimParams.conditions(ii,1) * stimParams.NIdaqRate)'),1,stimParams.numOfStimulators);
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
        % add fake tactile channel
        stimParams.vibrotactileStimulus = [stimParams.vibrotactileStimulus, ...
            zeros(size(stimParams.vibrotactileStimulus,1),1)];
        % Loop through the stimulus sequence
        for ii = 1:length(stimParams.stimOnsetsSamples)
            % Insert the tactile signal at the respective time points and stimulator
            stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii) + 1):int64(stimParams.stimOnsetsSamples(ii) + 1) +stimParams.triggerSamples - 1, ...
                stimParams.numOfStimulators + 1) = ones(stimParams.triggerSamples,1);
        end
end

%Save stimulus matrix, tsv info and make figures
fname = sprintf('%s_%s_%d.mat', stimParams.site, condition, runNum);

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

% Set TSV file information
onsets          = round(stimParams.stimOnsetsSecs,3);
duration        = stimParams.conditions(:,1);
trial_type      = stimParams.conditions(:,1);
trial_name      = stimParams.conditions(:,1);
stim_file       = repmat(fname, length(stimParams.stimOnsetsSecs),1);
stim_order      = stimParams.conditions(:,1);
stim_file_index = repmat('n/a', length(stimParams.stimOnsetsSecs),1);

stimParams.tsv = table(onsets, duration, trial_type, trial_name, stim_file, stim_file_index, stim_order);

% store for main script calls
stimulus.vibrotactileStimulus = stimParams.vibrotactileStimulus;

% Save
if ~exist(fullfile(vistadispRootPath, 'StimFiles',  fname), 'file')
    fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
    save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')
    save(fullfile(vistadispRootPath, 'StimFiles',  [sprintf('%s_%s', stimParams.site, condition),'_stimParams.mat']), 'stimParams');
else
    error('stimulus files for this experiment already exist, move or delete old files')
end

