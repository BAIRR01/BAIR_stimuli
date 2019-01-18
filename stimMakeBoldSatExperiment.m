function stimMakeBoldSatExperiment(stimParams, runNum, stimDurationSeconds, TR, movementRate, movementRateIndex)
% stimMakeTaskExperiment(stimParams,  runNum, TR)

site = stimParams.experimentSpecs.sites{1};

%For 0.33Hz, each trial consisted of three visual cues (the red circle
%turning green similar as in the first part, for 400 ms), presented at 3 s
%intercue interval, instructing the subject to close and open their hand,
%thus resulting in 0.33 Hz movement. In the second, third, and fourth
%tasks, respectively, 4, 7, or 13 cues were given with 2, 1, and 0.5 s
%intercue interval. Each trial was followed by a rest period until a total
%trial duration of 19 s was reached. Each task consisted of 13?15 trials of
%19 s.

frameRate             = stimParams.display.frameRate;
numberOfEventsPerRun  = 15; 
movementOnDuration    = 6;

% determine how many movements per trial for this movementRate
movementISI = 1/movementRate;
movementsPerTrial = round(movementOnDuration/movementISI);
if movementRateIndex == 1
    movementsPerTrial = movementsPerTrial+1;
end

% Determine trial length
switch(lower(stimParams.modality))
    case 'fmri'        
        %totalTrialDuration = 21; 
        movementOffDuration = 15;
        preScanPeriod       = 12;
    case {'ecog' 'eeg' 'meg'}
        %totalTrialDuration = 12; 
        movementOffDuration = 6;
        preScanPeriod       = 3;
end
preScanPeriod       = round(preScanPeriod/TR)*TR;
postScanPeriod      = preScanPeriod; % seconds
movementOffDuration = round(movementOffDuration/TR)*TR;
trialDuration = (stimDurationSeconds*movementsPerTrial)+(movementISI*(movementsPerTrial-1))+movementOffDuration;
%trialDuration  = round(totalTrialDuration/TR)*TR;

% One trial consists of multiple repeated movements
% Each movements has a stimDurationSeconds

% Generate onsets
onsets = cumsum([preScanPeriod ones([1 numberOfEventsPerRun-1])*trialDuration]);

% Round to multiples of TR
onsets = round(onsets/TR)*TR;

% Derive indices into the stimulus sequence (defined at frameRate)
onsetIndices  = round(onsets*(frameRate*.5))+1;

% Define total length of experiment
experimentLength = onsets(numberOfEventsPerRun)+trialDuration+postScanPeriod;

% Create a stimulus struct
stimulus = [];
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;

% No images are shown in this experiment, just fixation dot. So we just use
% the empty placeholder images from stimInitialize
stimulus.images     = stimParams.stimulus.images(:,:,end);

% Specify stimulus sequence at frame rate resolution
stimulus.seqtiming  = 0:(1/frameRate*2):experimentLength;
stimulus.seq        = ones(size(stimulus.seqtiming));

% For this motor task, the fixation changes are the task instructions
% (clench vs rest). So the fixation sequence is the one that needs to be
% filled with the movement cue onsets
stimulus.fixSeq     = stimulus.seq;

% Add the stimulus indices
% One trial consists of multiple repeated movements
% Each movements has a stimDurationSeconds

%movementISI = 1/movementRate;
%movementsPerTrial = round(movementOnDuration/movementISI)+1;
movementOnsets = 0:movementISI+stimDurationSeconds:(movementISI+stimDurationSeconds)*movementsPerTrial;

% Round to multiples of frameRate
movementOnsets  = round(movementOnsets*frameRate)/frameRate;
movementOnsetIndices  = round(movementOnsets*(frameRate*.5))+1;

imagesPerTrial = round((movementOnsets(end)+stimDurationSeconds)*(frameRate*.5));
imagesPerMovement = round(stimDurationSeconds*(frameRate*.5));

% Generate sequence of movements onsets per trial
sequencePerTrial = ones(1,imagesPerTrial);
sequencePerMovement = ones(1,imagesPerMovement);
for ii = 1:movementsPerTrial
    indices = movementOnsetIndices(ii) + (0:imagesPerMovement-1);
    sequencePerTrial(indices) = sequencePerMovement+1;
end

% Add all trials together
for ii = 1:numberOfEventsPerRun
    indices = onsetIndices(ii) + (0:length(sequencePerTrial)-1);
    stimulus.fixSeq(indices) = sequencePerTrial;
end

stimulus.onsets = onsets;

% Describe stimuli 
stimulus.cat            = [1 movementRateIndex+2]; %3,4,5 or 6
stimulus.categories     = {'REST', 'CLENCH'};
stimulus.movementRate   = movementRate; 

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri' 
        % No trigger sequence needed
    otherwise
         % Write binary trigger sequence:
        stimulus.trigSeq = zeros(size(stimulus.seqtiming));
        stimulus.trigSeq(onsetIndices) = stimulus.cat(2);
        stimulus.trigSeq(1)      = 255; %experiment onset
        stimulus.trigSeq(end)    = 255; %experiment offset
end

% Sparsification is OFF for motor because we need to sample data glove at
% rate very close to stimulus frame rate 
%maxUpdateInterval = 0.25;
%stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);

% Create stim_file name
fname = sprintf('%s_boldsat%d_%d.mat', site, movementRateIndex, runNum);

onset       = reshape(round(stimulus.onsets,3), [length(stimulus.onsets) 1]);
duration    = ones(length(stimulus.onsets),1) * stimDurationSeconds;
trial_type  = ones(length(stimulus.onsets),1);
trial_name  = repmat(stimulus.categories{2}, length(stimulus.onsets),1);
stim_file   = repmat(fname, length(stimulus.onsets),1);
stim_file_index = repmat('n/a', length(stimulus.onsets),1);
   
stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% Save
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end

