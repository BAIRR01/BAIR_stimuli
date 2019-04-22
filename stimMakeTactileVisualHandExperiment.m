function stimMakeTactileVisualHandExperiment(stimParams, runNum, directions, makeFigure)
% Makes a simple tactile hand experiment where the stimulators vibrate in
% several patterns and assumes one stimulator per finger.
%
%  Options are:
%
%     1. Ascending: stimulator 1 to 4, hopefully starting at index finger
%

% set some stimulus properties
stimulus = [];
stimulus.NIdaqRate              = stimParams.NIdaqRate;
stimulus.NIdaqNames             = stimParams.NIdaqNames;
stimulus.numOfStimulators       = stimParams.numOfStimulators;
stimulus.frameRate              = stimParams.display.frameRate;
stimulus.cmap                   = stimParams.stimulus.cmap;
stimulus.srcRect                = stimParams.stimulus.srcRect;
stimulus.dstRect                = stimParams.stimulus.destRect;
stimulus.display                = stimParams.display;

stimParams.fingerIdx              = [1:stimulus.numOfStimulators, cumsum(1:stimulus.numOfStimulators)];
stimParams.fingers                = {'index', 'middle', 'ring', 'little','all'};
stimParams.stimDurSamples         = round(stimParams.stimDurSecs * stimParams.NIdaqRate);
stimParams.stimDurFrames          = stimParams.stimDurSecs * stimParams.display.frameRate;

%condition based settings
if contains(lower(stimParams.condition),'weak')
    stimParams.tactileIntensity       = 1.25; % maximal value is 2, check state of amplifier, too
elseif contains(lower(stimParams.condition),'strong')
    stimParams.tactileIntensity       = 1.75; % maximal value is 2, check state of amplifier, too
else
    error('tactile strength undefined');
end

if contains(lower(stimParams.condition),'visualoff')
    stimParams.visualPulse            = false; % add a pulsating disk (true or false)
else
    stimParams.visualPulse            = true; % add a pulsating disk (true or false)
end

% warn if samples don't match secs
if (stimParams.stimDurFrames / stimParams.display.frameRate ~= stimParams.stimDurSecs)
    error('stimulus duration cannot be calculated in samples');
end

%% stimulus sequence timing
for jj = 1:length(directions)
    switch directions{jj}
        case {'Ascending','Descending'}
            %duration of the whole experiment
            stimParams.expDurSecs = stimParams.numSweeps * ... %number of sweeps
                (stimParams.stimDurSecs*stimParams.numOfStimulators ...%all fingers are stimulated during one sweep
                + stimParams.interSweepIntervalSecs)...%after each sweep there is a pause
                + stimParams.preScanDurSecs + stimParams.postScanDurSecs;%pauses at the beginning and the end of the experiment
            %durations of different (non-)stimulation periods
            stimParams.periodDursSecsVect  = [stimParams.preScanDurSecs, ... %pre-stimulation
                repmat(...
                [repmat(stimParams.stimDurSecs, 1, stimParams.numOfStimulators)... %one sweep of tactile stimulations
                stimParams.interSweepIntervalSecs],... % pause after each sweep
                1,stimParams.numSweeps)... % repeat for each sweep
                stimParams.postScanDurSecs];% add post-experimental pause
            %onsets of tactile(
            stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs, ... %first stimulus starts after pre-experiment pause
                repmat(...
                repmat(stimParams.stimDurSecs, 1, stimParams.numOfStimulators)... %during sweeps, one stimulus follows immediately after the other
                + [zeros(1,stimParams.numOfStimulators - 1) stimParams.interSweepIntervalSecs],... % add pause after each sweep to last duration
                1,stimParams.numSweeps-1)... % repeat for all but last sweep
                repmat(stimParams.stimDurSecs, 1, stimParams.numOfStimulators-1)]);%no stimulus starts after last stimuus
        case 'Blocked'
            %duration of the whole experiment
            stimParams.expDurSecs = stimParams.numSweeps * ... %number of sweeps
                (stimParams.stimDurSecs ...%one stimulus (all fingers) per sweep
                + stimParams.interSweepIntervalSecs)...%after each sweep there is a pause
                + stimParams.preScanDurSecs + stimParams.postScanDurSecs;%pauses at the beginning and the end of the experiment
            %durations of different (non-)stimulation periods
            stimParams.periodDursSecsVect  = [stimParams.preScanDurSecs, ... %pre-stimulation
                repmat(...
                [stimParams.stimDurSecs,... %one tactile stimulation
                stimParams.interSweepIntervalSecs],... % pause after each sweep of one stimulation
                1,stimParams.numSweeps)... % repeat for each sweep
                stimParams.postScanDurSecs];% add post-experimental pause
            %onsets of tactile(
            stimParams.stimOnsetsSecs = cumsum([stimParams.preScanDurSecs, ... %first stimulus starts after pre-experiment pause
                repmat(...
                stimParams.stimDurSecs... %one sweep = one stimulation
                + stimParams.interSweepIntervalSecs,... % add pause after each sweep = stimulation
                1,stimParams.numSweeps-1)... % repeat for all sweeps
                ]);
    end
end
% calculate in frames
stimParams.expDurFrames = round(stimParams.expDurSecs * stimParams.display.frameRate);%duration of the whole experiment in frames


stimParams.stimOnsetsSamples         = stimParams.stimOnsetsSecs*stimParams.NIdaqRate;%onsets of tactile stim in samples
stimParams.stimOnsetsSecsFrameNormed = round(stimParams.stimOnsetsSecs*stimParams.display.frameRate)/stimParams.display.frameRate;
stimParams.stimOnsetsFrames          = stimParams.stimOnsetsSecsFrameNormed*stimParams.display.frameRate;%onsets of new tactile stim in frames
stimParams.seqtiming                 = 0:(1/stimParams.display.frameRate):stimParams.expDurSecs;%frame onsets
stimParams.fixSeq                    = 55*ones(size(stimParams.seqtiming)); %determines whether resp. which type of fixation is shown
stimParams.seq                       = ones(size(stimParams.seqtiming));

stimulus.seqtiming = stimParams.seqtiming;
stimulus.fixSeq = stimParams.fixSeq;
stimulus.seq = stimParams.seq;

%% Make stimulus for the experiment

% signal for one tactile stimulus (e.g., per finger)
% carrier signal for one tactile stimulus
stimulusSignalBase = 0.5 + 0.5 * sin(-pi/2+linspace(0, ...
    (2*pi*stimParams.carrierFreq/stimParams.NIdaqRate * stimParams.stimDurSamples), ...
    stimParams.stimDurSamples)');
% modulating signal for one tactile stimulus
stimulusSignalModulator = stimParams.tactileIntensity/2 + stimParams.tactileIntensity/2 * sin(-pi/2+linspace(0, ...
    (2*pi*stimParams.modulatingFreq/stimParams.NIdaqRate * stimParams.stimDurSamples), ...
    stimParams.stimDurSamples)');
% modulating signal for visual stimulus
stimulusSignalModulatorFrames = 2 + sin(-pi/2+linspace(0, ...
    (2*pi*stimParams.modulatingFreq/stimParams.display.frameRate * stimParams.stimDurFrames), ...
    stimParams.stimDurFrames)');
% final signal for one tactile stimulus
stimParams.stimulusSignal = stimulusSignalBase .* stimulusSignalModulator;
% replace last entry with 0 to inactivate tactile stimulator
stimParams.stimulusSignal(end) = 0;
%initialize matrix to hold activation of each stimulator during each
%sample of the nidaq for one whole experiment
stimParams.vibrotactileStimulus = zeros(round(stimParams.expDurSecs*stimParams.NIdaqRate), stimParams.numOfStimulators);

%initialize matrix to hold position (x,y) and size of circle for each frame during experiment
stimParams.visualStimulus = zeros(round(stimParams.expDurSecs*stimParams.display.frameRate), 3);


%% prepare images
% X,Y center location for fixation
stimParams.fixX = (stimParams.stimulus.destRect(1)+stimParams.stimulus.destRect(3))/2;
stimParams.fixY = (stimParams.stimulus.destRect(2)+stimParams.stimulus.destRect(4))/2;
% store settings for hand image, used as fixation
stimParams.handImageDims = [0 0 900 900];
% store path to hand stimulus
stimParams.handImagePath = [fullfile(BAIRRootPath , 'handStimuliResources') '/Hand.png'];
% load the circle image
[circleImage,~,alpha] = imread([fullfile(BAIRRootPath , 'handStimuliResources') '/Circle.png']);
% add alpha values to the image
circleImage(:, :, 4) = alpha;
%initialize stimulus images to ensure dimension fit later calls
stimParams.images          = zeros([size(circleImage) 1]);
%add circle image
stimParams.images(:,:,:,1) = circleImage;

% store for use in main script
stimulus.handImageDims = stimParams.handImageDims;
stimulus.handImagePath = [fullfile(BAIRRootPath , 'handStimuliResources') '/Hand.png'];
stimulus.images = stimParams.images;

% check for visual congruency
if contains(lower(stimParams.condition),'visualcongruent')
    %set coordinates for each finger (index, middle, ring, little)
    stimParams.fingerCoords = [-8,-25; -23,-13; -27,4; -21, 21] * max(stimParams.handImageDims) * 0.01;
elseif contains(lower(stimParams.condition),'visualincongruent')
    stimParams.fingerCoords = [17,-22;-8,-25; -23,-13; -27,4] * max(stimParams.handImageDims) * 0.01;
elseif ~stimParams.visualPulse
    %set coordinates but tut the size of the circle will be 0
    stimParams.fingerCoords = [-8,-25; -23,-13; -27,4; -21, 21] * max(stimParams.handImageDims) * 0.01;
else
    error('visual condition undefined')
end

% check whether pulsing circle should be added
if stimParams.visualPulse
    %set base rectangle for circle
    circleBaseRect = 0.04 * CenterRectOnPointd(stimParams.handImageDims, stimParams.fixX, stimParams.fixY);
else
    %make the rectangle zero size
    %set base rectangle for circle
    circleBaseRect = 0.00 * CenterRectOnPointd(stimParams.handImageDims, stimParams.fixX, stimParams.fixY);
end
%pre-allocate array for circle coordinates
stimParams.dstRect2 = zeros(stimParams.expDurFrames + 1,4);

for jj = 1:length(directions)
    switch directions{jj}
        case 'Ascending'
            stimulatorOrder = repmat(1:stimParams.numOfStimulators, 1, stimParams.numSweeps)';
            stimulatorOrder2 = stimulatorOrder;
        case 'Descending'
            stimulatorOrder = repmat(stimParams.numOfStimulators:-1:1, 1, stimParams.numSweeps)';
            stimulatorOrder2 = stimulatorOrder;
        case 'Blocked'
            stimulatorOrder = repmat(1:stimParams.numOfStimulators, stimParams.numSweeps,1);
            stimulatorOrder2 = repmat(5, 1, stimParams.numSweeps);
    end
    
    % Loop through the stimulus sequence
    for ii = 1:length(stimParams.stimOnsetsSamples)
        % Insert the tactile signal at the respective time points and stimulator
        stimParams.vibrotactileStimulus(int64(stimParams.stimOnsetsSamples(ii)):int64(stimParams.stimOnsetsSamples(ii)) + stimParams.stimDurSamples - 1, ...
            stimulatorOrder(ii,:)) = repmat(stimParams.stimulusSignal,1,size(stimulatorOrder,2));
    end
    
    % Loop through the stimulus sequence
    for ii = 1:length(stimParams.stimOnsetsFrames)
        stimParams.visualStimulus(uint64(stimParams.stimOnsetsFrames(ii)):uint64(stimParams.stimOnsetsFrames(ii)) + stimParams.stimDurFrames - 1, 1) = stimParams.fingerCoords(stimulatorOrder(ii),1);
        stimParams.visualStimulus(uint64(stimParams.stimOnsetsFrames(ii)):uint64(stimParams.stimOnsetsFrames(ii)) + stimParams.stimDurFrames - 1, 2) = stimParams.fingerCoords(stimulatorOrder(ii),2);
        % Insert the tactile signal at the respective time points and stimulator
        stimParams.visualStimulus(uint64(stimParams.stimOnsetsFrames(ii)):uint64(stimParams.stimOnsetsFrames(ii)) + stimParams.stimDurFrames - 1, 3) = stimulusSignalModulatorFrames;
    end
    
    % Loop through every frame and set circle
    for ii = 1:stimParams.expDurFrames
        % position the rectangle on the finger
        stimParams.dstRect2(ii,:) = CenterRectOnPointd(circleBaseRect*stimParams.visualStimulus(ii,3), stimParams.fixX + stimParams.visualStimulus(ii,1), stimParams.fixY + stimParams.visualStimulus(ii,2));
    end
    
    switch lower(stimParams.modality)
        case 'fmri'
            % no trigger sequence needed
        otherwise
            % Write trigger sequence
    end
    
    %Save stimulus matrix, tsv info and make figures
    fname = sprintf('%s_%s%s_%d.mat', stimParams.site,stimParams.condition, directions{jj}(1:4), runNum);
    
    % check whether figures should be made
    if makeFigure
        % save figure with images of stimulus
        f = figure('visible', 'off');
        imagesc(stimParams.vibrotactileStimulus)
        title (sprintf('%s-%s', stimParams.condition, directions{jj}))
        xlabel('Stimulators');
        ylabel('Time (s)')
        xticks(stimParams.fingerIdx(1:4))
        xticklabels(stimParams.fingers(1:4))
        yticks(0:stimParams.NIdaqRate * 10:stimParams.expDurSecs * stimParams.NIdaqRate)
        yticklabels(0:10:stimParams.expDurSecs)
        saveas(f, fullfile(vistadispRootPath, 'StimFiles', sprintf('%s.png',fname(1:end-6))))
    end
    
    % Set TSV file information
    onsets          = round(stimParams.stimOnsetsSecs,3)';
    duration        = repmat(stimParams.stimDurSecs, length(stimParams.stimOnsetsSecs),1);
    trial_type      = stimParams.fingerIdx(stimulatorOrder2)';
    trial_name      = stimParams.fingers(stimulatorOrder2)';
    stim_file       = repmat(fname, length(stimParams.stimOnsetsSecs),1);
    stim_order      = repmat(directions(jj), length(stimParams.stimOnsetsSecs),1);
    stim_condition  = repmat(stimParams.condition, length(stimParams.stimOnsetsSecs),1);
    stim_file_index = repmat('n/a', length(stimParams.stimOnsetsSecs),1);
    
    stimParams.tsv = table(onsets, duration, trial_type, trial_name, stim_file, stim_file_index, stim_order, stim_condition);
    
    % store for main script calls
    stimulus.vibrotactileStimulus = stimParams.vibrotactileStimulus;
    stimulus.dstRect2 = stimParams.dstRect2;
    
    % Save
    if ~exist(fullfile(vistadispRootPath, 'StimFiles',  fname), 'file')
        fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
        save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')
        save(fullfile(vistadispRootPath, 'StimFiles',  [sprintf('%s_%s%s', stimParams.site,stimParams.condition, directions{jj}(1:4)),'_stimParams.mat']), 'stimParams');
    else
        error('stimulus files for this experiment already exist, move or delete old files')
    end
end

