function stimMakePRFExperiment(stimParams, runNum, stimulusDuration, isi, TR)

frameRate = stimParams.display.frameRate;

% Load aperture images
barApertures = stimMakePRFApertures(stimParams);

% Load carrier images
barCarriers   = stimMakePRFCarriers(stimParams);

% How many images shown for each aperture?
imagesPerPosition = round(TR / (stimulusDuration + isi));

% Total number of images
numberOfImages = imagesPerPosition * size(barApertures,3);

% Total time
totalTime = numberOfImages * (stimulusDuration + isi);

% Find a blank from the image sequence. This will be reused later.
flattenedApertures = reshape(barApertures, [], size(barApertures,3));
isblank = sum(flattenedApertures,1)==0;
blankImageIndex  = find(isblank,1);
 
% Image size
imageSizeInPixels = size(stimParams.stimulus.images);

% Recast to double centered at background color
barCarriers         = double(barCarriers);
backgroundColor     = mode(barCarriers(:));
barCarriers         = barCarriers  - backgroundColor;

% Make the images by multistimplying the carriers and the apertures
images = zeros(imageSizeInPixels(1), imageSizeInPixels(2), numberOfImages, 'double');
for ii = 1:numberOfImages
    idx = ceil(ii/imagesPerPosition);
    idx2 = randsample(size(barCarriers,3),1);
    images(:,:,ii) = barApertures(:,:,idx) .*  barCarriers(:,:,idx2);
end

% Recast as 8 bit integers
images = uint8(images + backgroundColor);

% Create a stimulus struct
stimulus = [];
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;
stimulus.images     = images;

% Create category labels for the apertures (will be used for trigger codes)
catNumbers          = ones(1, numberOfImages)*255; % BLANK 
catNumbers(1:28)    = 3:1:30;  % VERTICAL BARS   (moving left-right)
catNumbers(29:40)   = 59:1:70; % DIAGONAL BARS   (moving rightdown-leftup)
catNumbers(57:84)   = 31:1:58; % HORIZONTAL BARS (moving up-down)
catNumbers(85:96)   = 71:1:82; % DIAGONAL BARS   (moving leftdown-rightup)
catNumbers(113:140) = 30:-1:3; % VERTICAL BARS   (moving right-left)
catNumbers(141:152) = 83:1:94; % DIAGONAL BARS   (moving leftup-rightup)
catNumbers(169:196) = 58:-1:31;% HORIZONTAL BARS (moving down-up)
catNumbers(197:208) = 95:1:106;% DIAGONAL BARS   (moving rightup-leftdown)

% Create corresponding category labels for trial_name column in tsv file
% UGLY loops but didn't know how to do this otherwise...

catNames = cell(1,length(catNumbers));

stimId = 1;
for cc = 3:30
    catNames{cc} = ['VERTICAL-L-R-' num2str(stimId)]; stimId = stimId + 1;
end
stimId = 1;
for cc = 31:58
    catNames{cc} = ['HORIZONTAL-U-D-' num2str(stimId)]; stimId = stimId + 1;
end
stimId = 1;
for cc = 59:70
    catNames{cc} = ['DIAGONAL-RD-LU-' num2str(stimId)]; stimId = stimId + 1;
end
stimId = 1;
for cc = 71:82
    catNames{cc} = ['DIAGONAL-LD-RU-' num2str(stimId)]; stimId = stimId + 1;
end
stimId = 1;
for cc = 83:94
    catNames{cc} = ['DIAGONAL-LU-RD-' num2str(stimId)]; stimId = stimId + 1;
end
stimId = 1;
for cc = 95:106
    catNames{cc} = ['DIAGONAL-RU-LD-' num2str(stimId)]; stimId = stimId + 1;
end
catNames{255}    = 'BLANK';

% Add category labels to struct
stimulus.cat        = catNumbers;
stimulus.categories = catNames;

% Add duration and trialIndex into images
stimulus.duration   = ones(1,numberOfImages) * stimulusDuration;
stimulus.ISI        = ones(1,numberOfImages) * isi;
stimulus.trialindex = 1:1:numberOfImages;

% Create trial sequence
stimulus.seqtiming  = 0:1/frameRate:totalTime;
stimulus.seq        = zeros(size(stimulus.seqtiming))+blankImageIndex;

numIms    = round(stimulusDuration*frameRate);
numBlanks = round(isi*frameRate);
chunkSize = numIms + numBlanks;

for ii = 1:numberOfImages
    idx = (1:chunkSize) + (ii-1)*chunkSize;
    stimOrder = [repmat(ii, [1,numIms]) repmat(blankImageIndex, [1, numBlanks])];
    stimulus.seq(idx) = stimOrder;
    stimulus.onsets(ii) = stimulus.seqtiming(idx(1));
end

% Add pre and post stim period
switch(lower(stimParams.modality))
    case 'fmri'
        stimulus.prescan  = round(12/TR)*TR; % seconds
        stimulus.postscan = prescan; % seconds
    case {'ecog' 'eeg' 'meg'}      
        stimulus.prescan  = 3; % seconds
        stimulus.postscan = 3; % seconds 
    otherwise
        error('Unknown modality')
end

updatedSeqTimes = stimulus.seqtiming+stimulus.prescan;
prescanSeqTimes = 0:1/frameRate:stimulus.prescan-(1/frameRate);
postscanSeqTimes = (1/frameRate:1/frameRate:stimulus.postscan)+totalTime+stimulus.prescan;

stimulus.seqtiming = [prescanSeqTimes updatedSeqTimes postscanSeqTimes];
stimulus.seq = [zeros(1,length(prescanSeqTimes)) stimulus.seq zeros(1,length(postscanSeqTimes))];
stimulus.onsets = stimulus.onsets + stimulus.prescan;

% Get onsets
%[vals, Ai]          = unique(stimulus.seq);
%stimulus.onsets     = stimulus.seqtiming(Ai(vals~=0));

% Add fixation sequence
minDurationInSeconds = 1;
maxDurationInSeconds = 5;
fixSeq = createFixationSequence(stimulus, 1/frameRate, minDurationInSeconds, maxDurationInSeconds);
stimulus.fixSeq = fixSeq;

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri'
        % no triggers for fMRI
    otherwise
        % create an empty trigger sequence
        trigSeq  = zeros(size(stimulus.seq));
        % find the onsets of the stimuli in the sequence
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

% Create stimulus.mat filename
fname = sprintf('%s_prf_%d.mat', stimParams.site, runNum);

% Add table with elements to write to tsv file for BIDS
onset       = stimulus.onsets';
duration    = stimulus.duration';
ISI         = stimulus.ISI';
trial_type  = stimulus.cat(stimulus.trialindex)'; 
trial_name  = stimulus.categories(trial_type)';
stim_file   = repmat(fname, numberOfImages,1);
stim_file_index = stimulus.trialindex';

stimulus.tsv = table(onset, duration, ISI, trial_type, trial_name, stim_file, stim_file_index);

% remove empty cells from catNames
stimulus.categories = stimulus.categories(~cellfun('isempty',stimulus.categories));  

% save
stimulus.display  = stimParams.display;
stimulus.modality = stimParams.modality;
stimulus.site     = stimParams.site;

save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

