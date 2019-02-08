function stimMakeTactileHandExperiment(stimParams,  runNum, stimDuration, TR, numOfStimulators, directions , makeFigure)
% Makes a simple tactile hand experiment where the stimulators vibrate in
% several patterns and assumes one stimulator per finger.
% 
%  Options are:
%
%     1. Ascending: Thumb to little finger (depending on setup)
%     2. Descending: Little finger to thumb
%     3. Random: A random sequence of the stimulators repeated each cycle
%     4. Staircase1: Thumb to little finger followed by little finger to thumb
%     5. Staircase2: Little finger to thumb followed by thumb to little finger
%     
%     Note: Each case will have no blank/break in stimulation between
%     stimulator, but will have a blank/break (determined by the ISI
%     given)between sweeps
%
    
switch(lower(stimParams.modality))
    case 'fmri'
        nrCycles            = 6; % how many sweeps
        vibrFreq            = 30;
        preScanPeriod       = round(12/TR)*TR; % seconds
        pulseOnDur          = 400; %msec
        nrPulsesPerStim     = 12; % within each stimulus duration
        stimDuration        = round(stimDuration/TR)*TR;
        isi                 = stimDuration;
        
    case {'ecog' 'eeg' 'meg'}
        isi                 = stimDuration;
        nrCycles            = 6; % how many sweeps
        vibrFreq            = 30;
        preScanPeriod       = 6; % seconds
        pulseOnDur          = 400; % msec
        nrPulsesPerStim     = 12; % within each stimulus duration
    otherwise
        error('Unknown modality')
end

% Find the number of events we can fit in desired experiment time and
postScanPeriod   = preScanPeriod;
experimentLength = nrCycles*(stimDuration*numOfStimulators + isi) + preScanPeriod + postScanPeriod;
allISIs          = repmat([zeros(1,numOfStimulators-1) isi],1,nrCycles);
frameRate        = stimParams.display.frameRate;
onsets           = cumsum([preScanPeriod stimDuration + allISIs(1:end-1)]);
onsets           = round(onsets*frameRate)/frameRate;
onsetFrameIdx    = round(onsets*(frameRate*.5));
onsetSampleIdx   = onsets*stimParams.NIdaqRate;

% set some stimulus properties
stimulus = [];
stimulus.cat        = [1 2 3 4 5];
stimulus.categories = {'thumb', 'index', 'middle', 'ring', 'little'};
stimulus.prescan    = preScanPeriod;
stimulus.postscan   = postScanPeriod;
stimulus.onsets     = onsets;
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;
stimulus.vibrFreq   = vibrFreq;
stimulus.NIdaqRate  = stimParams.NIdaqRate;
stimulus.NIdaqNames = stimParams.NIdaqNames;
stimulus.numCycles  = nrCycles;

stimulus.seqtiming  = 0:(1/frameRate)*2:experimentLength;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));
stimulus.seq        = ones(size(stimulus.seqtiming));

% make a blank to insert between simulus presentations
screenRect          = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));
images              = ones([screenRect 1], 'uint8');
images(:,:,1)       = 127;
stimulus.images     = images;

%% Figure out duration of On/Off periods and make the stimulus 

%figure out the time between pulses and convert to msec
pulseOffDur = round(((stimDuration*1000)/ nrPulsesPerStim) - pulseOnDur);

% Store durations in samples
pulseOnDurSamp = (pulseOnDur/1000) * stimParams.NIdaqRate;
pulseOffDurSamp = (pulseOffDur/1000) * stimParams.NIdaqRate;

% One base unit of tactile stimulation vibrating with tactile frequency for
% base on time and not vibrating for tactile off time (in samples)
vibrInSamples = vibrFreq/ stimParams.NIdaqRate;
pulse = [1 + sin(linspace(0, 2 * pi * vibrInSamples * pulseOnDurSamp ,...
    pulseOnDurSamp)'); zeros(pulseOffDurSamp, 1)];

% signal for one tactile stimulus (e.g., per finger)
stimulusSignal          = repmat(pulse, nrPulsesPerStim, 1); % signal for one tactile stimulus
stimulusSignalDuration  = length(stimulusSignal);%in samples

vibrotactileStimulus = zeros(round(experimentLength*stimParams.NIdaqRate), numOfStimulators);

for ii = 1:length(directions)
    switch directions{ii}
        case 'ascending'
            stimulatorOrder = 1:numOfStimulators;
            experimentOrder = repmat(stimulatorOrder, 1, nrCycles);
        case 'descending'
            stimulatorOrder = flip(1:numOfStimulators);
            experimentOrder = repmat(stimulatorOrder, 1, nrCycles);
        case 'random'
            stimulatorOrder =  randperm(numOfStimulators);
            experimentOrder = repmat(stimulatorOrder, 1, nrCycles);
        case {'staircase1' , 'staircase2'}
            if contains (directions{ii}, '1')
                stimulatorOrder = [1:numOfStimulators flip(1:numOfStimulators)];
            elseif contains (directions{ii}, '2')
                stimulatorOrder = [flip(1:numOfStimulators) 1:numOfStimulators];
            end
            if mod(nrCycles, 2) ~= 0
                experimentOrder = [repmat(stimulatorOrder, 1, floor(nrCycles/2)) stimulatorOrder(1:numOfStimulators)];
            else
                experimentOrder = repmat(stimulatorOrder, 1, nrCycles/2);
            end
    end
    
    % Loop through the stimulus sequence
    for jj = 1:length(onsets)
        % Insert the tactile signal at the respective time point and stimulator
        vibrotactileStimulus(onsetSampleIdx(jj):onsetSampleIdx(jj)+ stimulusSignalDuration -1, ...
             experimentOrder(jj))  = stimulusSignal;
    end
  stimulus.vibrotactileStimulus = vibrotactileStimulus;
  
    switch lower(stimParams.modality)
        case 'fmri'
            % no trigger sequence needed
        otherwise
            % Write trigger sequence
            stimulus.trigSeq        = zeros(length(stimulus.fixSeq),1);
            stimulus.trigSeq(onsetFrameIdx) =  stimulus.cat(experimentOrder);
            stimulus.trigSeq(1)     = 255; %experiment onset
            stimulus.trigSeq(end)   = 255; %experiment offset
    end
    
    %Save stimulus matrix, tsv info and make figures
    fname = sprintf('%s_%sHand_%d.mat', stimParams.site,directions{ii}, runNum);
    
    % check whether figures should be made
    if makeFigure == 1
        % save figure with images of stimulus
        f = figure('visible', 'off');
        imagesc(vibrotactileStimulus)
        title (sprintf('Hand-%s', directions{ii}))
        xlabel('Stimulators');
        ylabel(sprintf('Time (%d samples/s)', stimParams.NIdaqRate));
        saveas(f, fullfile(vistadispRootPath, 'StimFiles', sprintf('%s.png',fname)))
    end
    
    % Set TSV file information
    onsets          = round(stimulus.onsets,3)';
    duration        = repmat(stimDuration, length(onsets),1);
    trial_type      = stimulus.cat(experimentOrder)';
    trial_name      = stimulus.categories(experimentOrder)';
    stim_file       = repmat(fname, length(stimulus.onsets),1);
    stim_order      = repmat(directions(ii), length(onsets),1);
    stim_file_index = repmat('n/a', length(stimulus.onsets),1);
    
    stimulus.tsv = table(onsets, duration, trial_type, trial_name, stim_file, stim_file_index, stim_order);
    
    % Save
    fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
    save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')
end

