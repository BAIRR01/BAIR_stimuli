function stimMakeHRFExperiment(stimParams, runNum, stimDurationSeconds, onsetTimeMultiple, stimulusType, TR)
%stimMakeHRFExperiment(stimParams, runNum, stimDurationSeconds, onsetTimeMultiple, stimulusType)
%
%  Stimuli are presented for stimDurationSeconds, with an exponentially
%  distributed interstimulus interval (mean ~9s, range 3-24s).
% 
%  stimulusType = ... 
%     {'hrfPattern' 'hrfPatternInverted' 'hrfChecker' 'hrfCheckerInverted'};
    
% See s_Make_BAIR_Visual_Experiments for context


%% Make the images

% Determine if we're creating the master or loading resizing for a specific display
site = stimParams.experimentSpecs.Row{1};
frameRate = stimParams.display.frameRate;

switch site
    case 'Master'

        % Experiment specs      
        if strcmpi(stimulusType, 'HRFPATTERNBREATHINGCHALLENGE')
            numberOfEventsPerRun  = 71;
        else
            numberOfEventsPerRun  = 32;
        end
        
        preScanPeriod         = round(30/TR)*TR;
        postScanPeriod        = preScanPeriod;
        minimumISIinSeconds   = 3;
        maximumISIinSeconds   = 24;
       
        % Specify the pattern density in some arbitrary unit: a value of 20
        % gives a moderately dense pattern: See Kay et al, 2013 PLOS CB, figure 5B
        stripesPerImage = 20;
        imageSizeInPixels = size(stimParams.stimulus.images);
        
        % The number of images is 2 x the number of stimuli + 1. We multiply by 2
        % because the stimuli may be shown as paired images (a pattern and its
        % contrast reversal). We add 1 to define the blank stimulus.
        numImages = numberOfEventsPerRun*2+1;
        
        blankImageIndex = numImages;
        
        backgroundIntensity = 128;    
        images = zeros([imageSizeInPixels numImages], 'uint8')+backgroundIntensity;

        % Loop to make the images
        for ii = 1:numberOfEventsPerRun
            
            disp(['Creating stimulus number: ' num2str(ii)]);
            
            if contains(stimulusType, 'checker','IgnoreCase',true)
                imageForThisTrial = createCheckerboard(stimParams, round(stripesPerImage/3));

            elseif contains(stimulusType, 'pattern','IgnoreCase',true)
                imageForThisTrial = createPatternStimulus(stimParams, stripesPerImage);

            else 
                error('Unknown stimulus type %s', stimulusType);

            end

            % Double to unsigned 8 bit integer, needed for vistadisp
            image8Bit             = uint8((imageForThisTrial+.5)*255);

            % The contrast reversed images will be used for only those experiments
            % where the pulse has a contrast reversal, otherwise ignored
            imageContrastReversed = uint8((-imageForThisTrial+.5)*255);

            images(:,:,ii) = image8Bit;
            images(:,:,ii+numberOfEventsPerRun) = imageContrastReversed;
            
        end
        
        % Store the images in the stimulus structure used by vistadisp
        stimulus = [];
        stimulus.cmap         = stimParams.stimulus.cmap;
        stimulus.srcRect      = stimParams.stimulus.srcRect;
        stimulus.dstRect      = stimParams.stimulus.destRect;
        stimulus.display      = stimParams.display;
        stimulus.categories   = {upper(stimulusType)};
        stimulus.images       = images;
        
        % Add a category label to be able to combine stimuli across runs
        stimulus.cat          = ones(1,numberOfEventsPerRun);
        stimulus.duration     = ones(1,numberOfEventsPerRun)*stimDurationSeconds;
        
            % Specify event timing, with an exponential distribution of ISIs
        [onsets, onsetIndices] = getExponentialOnsets(numberOfEventsPerRun, preScanPeriod, ...
            minimumISIinSeconds, maximumISIinSeconds, onsetTimeMultiple, frameRate);
        stimulus.onsets       = onsets;
        
        % Define total length of stimulation sequence at frame rate resolution
        stimulus.seqtiming    = 0:1/frameRate:onsets(numberOfEventsPerRun)+stimDurationSeconds+postScanPeriod;
        stimulus.seq          = zeros(size(stimulus.seqtiming))+blankImageIndex;
        
        % Add fixation sequence
        minDurationInSeconds = 1;
        maxDurationInSeconds = 5;
        fixSeq = createFixationSequence(stimulus, 1/frameRate, minDurationInSeconds, maxDurationInSeconds);
        stimulus.fixSeq = fixSeq;
        
        % Randomize the assignment of image to trial
        trialIndex = randperm(numberOfEventsPerRun);
        stimulus.trialindex = trialIndex;
        stimulus.seq(onsetIndices) = trialIndex;
        
        % Specify the image sequence within each stimulus, which depends on the
        % stimulus duration and the dwell time per image, as well as whether or not
        % we include a contrast reversal
        imagesPerTrial = round(stimDurationSeconds*frameRate);
        sequencePerTrial = zeros(1,imagesPerTrial);

        switch lower(stimulusType)
            case {'hrfpatterninverted' 'hrfcheckerinverted'}
                % paired images per trial, ie an image and its contrast reversal
                contrastReversal = true;
            otherwise
                % single image pre trial
                contrastReversal = false;
        end

        % Add the contrast reversed stimuli to the sequence
        if contrastReversal
            sequencePerTrial(imagesPerTrial/2+1:imagesPerTrial) = numberOfEventsPerRun;
        end

        % Add the stimulus sequence index
        for ii = 1:numberOfEventsPerRun
            indices = onsetIndices(ii) + (0:imagesPerTrial-1);
            stimulus.seq(indices) = sequencePerTrial + trialIndex(ii);
        end
    
         % Generate a save name 
        fname = sprintf('%s_%s_%d.mat', site, lower(stimulusType), runNum);
        
    otherwise        
        % Resize the Master stimuli to the required stimulus size for this
                % modality and display
        fprintf('[%s]: Loading Master stimuli for: %s\n', mfilename, site);

        % Load the Master stimuli
        stimulus = loadBAIRStimulus(lower(stimulusType), 'Master', runNum);
        
        % The desired new ImageSize
        imageSizeInPixels = size(stimParams.stimulus.images);
        
        % Resize         
        fprintf('[%s]: Resizing Master stimuli for: %s\n', mfilename,  site);
        images = imresize(stimulus.images, imageSizeInPixels);
        
        % Soft circular mask (1 pixel of blurring per 250 pixels in the image)
        supportDiameter       = imageSizeInPixels(1);
        maskRadius            = (stimParams.stimulus.srcRect(3) - stimParams.stimulus.srcRect(1))/2;
        circularMask          = mkDisc(supportDiameter, maskRadius, (imageSizeInPixels+1)./2, 1/250 * imageSizeInPixels(1));
        imagesDouble          = double(images)/255-.5;
        imagesMasked          = bsxfun(@times,imagesDouble,  circularMask);
        stimulus.images       = uint8((imagesMasked+.5)*255);
        
        % Overwrite master stimulus with settings for this modality
        % and display
        stimulus.cmap         = stimParams.stimulus.cmap;
        stimulus.srcRect      = stimParams.stimulus.srcRect;
        stimulus.dstRect      = stimParams.stimulus.destRect;
        stimulus.display      = stimParams.display;

        % Add triggers for non-fMRI modalities
        switch lower(stimParams.modality)
            case 'fmri'
                % no trigger sequence needed
            otherwise
                
                % Create an empty trigger sequence
                trigSeq            = zeros(size(stimulus.seq));
                blankImageIndex    = mode(stimulus.seq);
                % Write image identity into trigger sequence:
                stimulus.trigSeq   = stimulus.seq;
                stimulus.trigSeq(stimulus.seq == blankImageIndex) = 0;
                % Find the onsets of the stimuli in the sequence
                [~,onsetIndices] = intersect(round(stimulus.seqtiming,4),round(stimulus.onsets,4));
                assert(length(onsetIndices) == length(stimulus.onsets));
                % use the CATEGORICAL labels as trigger codes
                trigSeq(onsetIndices) = stimulus.cat(stimulus.seq(onsetIndices));
                % add task ONSET and OFFSET trigger
                trigSeq(1)   = 256;
                trigSeq(end) = 256;
                stimulus.trigSeq = trigSeq;
        end
        
        % Sparsify
        maxUpdateInterval = 0.25;
        stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);
        
        % Generate a save name 
        fname = sprintf('%s_%s_%d.mat', site, lower(stimulusType), runNum);
        
        % Add table with elements to write to tsv file for BIDS
        onset       = reshape(round(stimulus.onsets,3), [length(stimulus.onsets) 1]);
        duration    = ones(length(stimulus.onsets),1) * stimDurationSeconds;
        trial_type  = ones(length(stimulus.onsets),1);
        trial_name  = repmat(upper(stimulusType), length(stimulus.onsets),1);
        stim_file   = repmat(fname, length(stimulus.onsets),1);
        stim_file_index = reshape(stimulus.trialindex, [length(stimulus.onsets) 1]);
        
        stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);
  
        % add modality
        stimulus.modality = stimParams.modality;    
       
end

% Save
stimulus.site     = site;

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
