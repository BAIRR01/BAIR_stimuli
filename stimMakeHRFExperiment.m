function stimMakeHRFExperiment(stimParams, runNum, stimDurationSeconds, dwellTimePerImage, stimulusType)
%stimMakeHRFExperiment(stimParams, runNum, stimDurationSeconds, dwellTimePerImage, stimulusType)
%
%  Stimuli are presented for stimDurationSeconds, with an exponentially
%  distributed interstimulus interval (mean ~9s, range 3-24s).
% 
%  stimulusType = ... 
%     {'hrfPattern' 'hrfPatternInverted' 'hrfChecker' 'hrfCheckerInverted'};
    
% See s_Make_BAIR_Visual_Experiments for context

%% Experiment timing

numberOfEventsPerRun  = 32;  
minimumISIinSeconds   = 3; 
maximumISIinSeconds   = 24;

% Specify event timing, with an exponential distribution of ISIs
[onsets, onsetIndices] = getExponentialOnsets(numberOfEventsPerRun,...
    minimumISIinSeconds, maximumISIinSeconds, dwellTimePerImage);

%% Make the images

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

% determine if we're creating the master or loading resizing for a specific display
site = stimParams.experimentSpecs.Row{1};

switch site
    case 'Master'

        % Loop to make the images
        for ii = 1:numberOfEventsPerRun
            
            disp(['Creating stimulus number: ' num2str(ii)]);
            
            if contains(stimulusType, 'checker','IgnoreCase',true)
                imageForThisTrial = createCheckerboard(stimParams, stripesPerImage*2);

            elseif contains(stimulusType, 'pattern','IgnoreCase',true)
                imageForThisTrial = createPatternStimulus(stimParams, stripesPerImage);

            else 
                error('Unknown stimulus type %s', stimulusType);

            end

            % Soft circular mask (1 pixel of blurring per 250 pixels in the image)
            circularMask          = mkDisc(imageSizeInPixels, imageSizeInPixels(1)/2, (imageSizeInPixels+1)./2, 1/250 * imageSizeInPixels(1));    
            imageForThisTrial     = imageForThisTrial .* circularMask;

            % Double to unsigned 8 bit integer, needed for vistadisp
            image8Bit             = uint8((imageForThisTrial+.5)*255);

            % The contrast reversed images will be used for only those experiments
            % where the pulse has a contrast reversal, otherwise ignored
            imageContrastReversed = uint8((-imageForThisTrial+.5)*255);

            images(:,:,ii) = image8Bit;

            images(:,:,ii+numberOfEventsPerRun) = imageContrastReversed;
        end
        
    otherwise        
        % Resize the Master stimuli to the required stimulus size for this
                % modality and display
        disp(['Loading Master stimuli for: ' site]);

        % Load the Master stimuli
        masterImages = loadBAIRStimulus(stimulusType, 'Master', runNum);
    
                 
        disp(['Resizing Master stimuli for: ' site]);
        images = imresize(masterImages, size(stimParams.stimulus.images));
end


% This is the stimulus structure used by vistadisp
stimulus = [];
stimulus.cmap         = stimParams.stimulus.cmap;
stimulus.srcRect      = stimParams.stimulus.srcRect;
stimulus.dstRect      = stimParams.stimulus.destRect;
stimulus.images       = images;
stimulus.seqtiming    = 0:dwellTimePerImage:300;
stimulus.seq          = zeros(size(stimulus.seqtiming))+blankImageIndex;

% Randomize the assignment of image to trial
imageIndex = randperm(numberOfEventsPerRun);
stimulus.seq(onsetIndices) = imageIndex;

% Specify the image sequence within each stimulus, which depends on the
% stimulus duration and the dwell time per image, as well as whether or not
% we include a contrast reversal
imagesPerTrial = round(stimDurationSeconds/dwellTimePerImage);
sequencePerTrial = zeros(1,imagesPerTrial);

switch lower(stimulusType)
    case {'pattern' 'checker'}
        % single image pre trial
        contrastReversal = false;
    case {'patterninverted' 'checkerinverted'}
         % paired images per trial, ie an image and its contrast reversal
        contrastReversal = true;
end

% Add the contrast reversed stimuli to the sequence
if contrastReversal
    sequencePerTrial(imagesPerTrial/2+1:imagesPerTrial) = numberOfEventsPerRun;
end

for ii = 1:numberOfEventsPerRun
    indices = onsetIndices(ii) + (0:imagesPerTrial-1);
    stimulus.seq(indices) = sequencePerTrial + imageIndex(ii);
end

% Create the fixation dot color change sequence
stimulus.fixSeq       = ones(size(stimulus.seqtiming));

this_frame = 0;
minDurationInSeconds = 1;
maxDurationInSeconds = 5;

minDurationInImageNumber = round(minDurationInSeconds / dwellTimePerImage);
maxDurationInImageNumber = round(maxDurationInSeconds / dwellTimePerImage);

while true
    % wait between minDurationInSeconds and maxDurationInSeconds before
    % flipping the dot color
    isi = randperm(maxDurationInImageNumber-minDurationInImageNumber,1)+minDurationInImageNumber-1;
    this_frame = this_frame + isi;
    if this_frame > length(stimulus.fixSeq), break; end
    stimulus.fixSeq(this_frame:end) = mod(stimulus.fixSeq(this_frame-1),2)+1;
end

% add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri' 
    otherwise
        stimulus.trigSeq  = double(stimulus.seq>0);
        stimulus.diodeSeq = stimulus.trigSeq;
end

% create stimulus.mat filename
fname = sprintf('hrf%s_%s_%d', stimulusType, site, runNum);

% add table with elements to write to tsv file for BIDS
onset = reshape(onsets, [numberOfEventsPerRun 1]);
duration = ones(numberOfEventsPerRun,1) * stimDurationSeconds;
trial_type = repmat('Pulse', numberOfEventsPerRun,1);
stim_file = repmat(fname, numberOfEventsPerRun,1);
stim_file_index = reshape(imageIndex, [numberOfEventsPerRun 1]);

stimulus.tsv = table(onset, duration, trial_type, stim_file, stim_file_index);

% NOTES:
% columns here were determined based on example script and file in:
% /Volumes/server/Projects/BAIR/MRI/data/visual/wl_subj001_20170823/code/makeHRFTSVFiles.m
% /Volumes/server/Projects/BAIR/BIDS_example/visualFullSet/sub-wlsubj014/ses-nyu3T01/func/sub-wlsubj014_ses-nyu3T01_task-spat_run-01_bold.tsv
% TO WRITE TSV FILE OUT use this function: 
% writetable(stimulus.tsv, 'output.tsv', 'FileType','text', 'Delimiter', '\t') 

% QUESTIONS: 
% Q1: should onsets be the ones derived with getExponentialOnsets (below)?
%     or should they be in units of dwellTimePerImage? In that case use this:
%           blankIdx = mode(stimulus.seq);
%           eventIdx = stimulus.seq ~= blankIdx;
%           onset = stimulus.seqtiming(eventIdx)
%           stim_file_index = stimulus.seq(eventIdx)
%     !!! onsets from getExponentialOnsets appears to be shifted relative to using
%         stimulus.seq? !!!
% Q3: are these the columns we want in the TSV file?
% Q4: writing out the TSV file should happen in stimulus presentation code? e.g. BAIR_FMRI
%     (probably needs to be updated anyway to accomodate filenaming changes?)
% Q5: do we want to add StimParams to stimulus file?

% save

stimulus.display  = stimParams.display;
stimulus.modality = stimParams.modality;
stimulus.site     = site;
stimulus          = bairCheckStimulus(stimulus);

save(fullfile(vistadispRootPath, 'Retinotopy', 'storedImagesMatrices',  fname), 'stimulus')

end

function [onsets, indices] = getExponentialOnsets(numStimuli, minISI, maxISI,temporalResolution)

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
ISIs = round(ISIs/temporalResolution)*temporalResolution;

% Compute the cummulative sum of ISIs to get the onset times
onsets = cumsum([8 ISIs(randperm(numStimuli-1))]);
indices  = round(onsets/temporalResolution);
% % Debug
% figure(2), clf; set(gcf, 'Color', 'w')
% set(gca, 'FontSize', 24, 'XTick', 0:60:300, 'YTick', []); hold on
% stem(onsets, ones(1,numStimuli+1), 'LineWidth', 2)
% xlabel('Time (s)')
% hgexport(gcf, fullfile(BAIRRootPath, 'figures', 'HRF_onsets.eps'))

end
