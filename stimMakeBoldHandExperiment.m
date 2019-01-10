function stimMakeBoldHandExperiment(stimParams, runNum, stimDurationSeconds, onsetTimeMultiple, TR)
% stimMakeTaskExperiment(stimParams,  runNum, TR)

site = stimParams.experimentSpecs.sites{1};

%  After the first 30 s rest, the circle turned green 54 times for 500 ms
%  at intertrial intervals varying from 3 to 18 s, with an average of 8 s.

numberOfEventsPerRun = 54;
frameRate            = stimParams.display.frameRate;

preScanPeriod         = round(30/TR)*TR;
postScanPeriod        = preScanPeriod;
minimumISIinSeconds   = 3;
maximumISIinSeconds   = 18;

% Generate onsets (same as for visual HRF experiment
[onsets, onsetIndices] = getExponentialOnsets(numberOfEventsPerRun, preScanPeriod, ...
    minimumISIinSeconds, maximumISIinSeconds, onsetTimeMultiple, frameRate);

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

% Describe stimuli
stimulus.cat        = [1 2];
stimulus.categories = {'REST', 'CLENCH'};

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri' 
        % no trigger sequence needed
    otherwise
        % Write binary trigger sequence:
        stimulus.trigSeq = stimulus.fixSeq-1; % no trigger should equal zero
end

% Sparsify stimulus sequence % ADAPT THIS FOR MOTOR??
%maxUpdateInterval = 0.25;
%stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);

% If we want to use different ISIs for ECOG, use switch like below:
% switch(lower(stimParams.modality))
%     case 'fmri'
%         minISI   = 3; % seconds
%         maxISI   = 18;  % seconds
%         prescan  = round(30/TR)*TR; % seconds 
%         postscan = prescan; 
% 
%         % Jitter ITIs
%         ISIs = linspace(ISI_min,ISI_max,numStim-1);                
% 
%         % Round off to onsetMultiple
%         %ISIs = round(ISIs/onsetTimeMultiple)*onsetTimeMultiple;
%         % Round off to TR
%         ISIs = round(ISIs/TR)*TR;
%         
%     case {'ecog' 'eeg' 'meg'}
%         ITI_min  = 1.25;
%         ITI_max  = 1.75;
%         prescan  = 3; % seconds
%         postscan = 3; % seconds
% 
%         % Jitter ITIs
%         ITIs = linspace(ITI_min,ITI_max,numberOfStimuli-1);
% 
%     otherwise
%         error('Unknown modality')
% end


% Create stim_file name
fname = sprintf('%s_boldhand_%d.mat', site, runNum);

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