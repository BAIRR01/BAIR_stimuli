function stimMakeBoldSatExperiment(stimParams, runNum, stimDurationSeconds, onsetTimeMultiple, TR, movementRate)
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
preScanPeriod         = round(30/TR)*TR;
postScanPeriod        = preScanPeriod;

% REPLACE THIS WITH INTERVALS BASED ON MOVEMENT RATE
%numberOfEventsPerRun  = 54;
%minimumISIinSeconds   = 3;
%maximumISIinSeconds   = 18;

% Generate onsets (same as for visual HRF experiment
%[onsets, onsetIndices] = getExponentialOnsets(numberOfEventsPerRun, preScanPeriod, ...
%    minimumISIinSeconds, maximumISIinSeconds, onsetTimeMultiple, frameRate);

numberOfEventsPerRun  = 15; 
eventLength           = 19; % seconds

% Define total length of experiment
experimentLength = onsets(numberOfEventsPerRun)+stimDurationSeconds+postScanPeriod;

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
stimulus.seqtiming  = 0:1/frameRate:experimentLength;
stimulus.seq        = ones(size(stimulus.seqtiming));

% For this motor task, the fixation changes are the task instructions
% (clench vs rest). So the fixation sequence is the one that needs to be
% filled with the trial onsets
stimulus.fixSeq     = stimulus.seq;

% Add the stimulus indices
imagesPerTrial = round(stimDurationSeconds*frameRate);
sequencePerTrial = ones(1,imagesPerTrial);

for ii = 1:numberOfEventsPerRun
    indices = onsetIndices(ii) + (0:imagesPerTrial-1);
    stimulus.fixSeq(indices) = sequencePerTrial+1;
end

stimulus.onsets = onsets;

% Describe stimuli UPDATE TO REFLECT MOVEMENTRATE SO DATA CAN BE
% CONCATENATED LATER (analogous to visual stimuli category coding)
stimulus.cat        = [1 2];
stimulus.categories = {'REST', 'CLENCH'};
stimulus.rateinHz   = movementRate;

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri' 
        % no trigger sequence needed
    otherwise
        % Write binary trigger sequence:
        stimulus.trigSeq = stimulus.fixSeq-1; % no trigger should equal zero
end

% Sparsification is OFF for motor because we need to sample data glove at
% rate very close to stimulus frame rate 
%maxUpdateInterval = 0.25;
%stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);

% Create stim_file name
fname = sprintf('%s_boldsat%dHz_%d.mat', movementRate, site, runNum);

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


function [onsets, indices] = getExponentialOnsets(numStimuli, preScanPeriod, minISI, maxISI, onsetTimeMultiple, frameRate)

% Draw numStimuli ISIs from an exponential distribution from [minISI maxISI]
x = linspace(0,1,numStimuli-1) * (1-exp(-minISI)) + exp(-minISI);
ISIs = -log(x)/minISI;
ISIs = ISIs*(maxISI-minISI)+minISI;

% % DEBUG
% disp(mean(ISIs))
%
% figure(1),clf, set(gcf, 'Color', 'w')
% set(gca, 'FontSize', 24); hold on
% plot(ISIs, 'o-', 'LineWidth', 4, 'MarkerSize', 12); axis tight
% ylabel('ITI (s)'); xlabel('Trial')

% Round off the ISIs to multiples of temporalResolution
ISIs = round(ISIs/onsetTimeMultiple)*onsetTimeMultiple;

% Compute the cumulative sum of ISIs to get the onset times
prescan  = round(preScanPeriod/onsetTimeMultiple)*onsetTimeMultiple;
onsets   = cumsum([prescan ISIs(randperm(numStimuli-1))]); 

% Match the stimulus presentation to the frame rate
onsets   = round(onsets*frameRate)/frameRate;

% Derive indices into the stimulus sequence (defined at temporalResolution)
indices  = round(onsets*frameRate)+1;

% % Debug
% figure(2), clf; set(gcf, 'Color', 'w')
% set(gca, 'FontSize', 24, 'XTick', 0:60:300, 'YTick', []); hold on
% stem(onsets, ones(1,numStimuli+1), 'LineWidth', 2)
% xlabel('Time (s)')
% hgexport(gcf, fullfile(BAIRRootPath, 'figures', 'HRF_onsets.eps'))

end