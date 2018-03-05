function stimMakePRFExperiment(stimParams, runNum, stimulusDuration, isi)

frameRate = stimParams.display.frameRate;

% Load aperture images
barApertures = stimMakePRFApertures(stimParams);

% Load carrier images
barCarriers   = stimMakePRFCarriers(stimParams);

% How many images shown for each aperture?
imagesPerPosition = round(1 / (stimulusDuration + isi));

% Total number of images
numberOfImages = imagesPerPosition * size(barApertures,3);

% Total time
totalTime = numberOfImages * (stimulusDuration + isi);

% Find a blank from the image sequence. This will be reused later.
flattenedApertures = reshape(barApertures, [], size(barApertures,3));
isblank = sum(flattenedApertures,1)==0;
blankImageIndex  = find(isblank,1);
% 
% % Total time points
% numberOfTimePoints = 1/dwellTimePerImage * size(barApertures,3);

% Image size
imageSizeInPixels = size(stimParams.stimulus.images);

% Recast to double centered at background color
barCarriers         = double(barCarriers);
backgroundColor     = mode(barCarriers(:));
barCarriers         = barCarriers  - backgroundColor;

% Make the images by mulitplying the carriers and the apertures
images = zeros(imageSizeInPixels(1), imageSizeInPixels(2), numberOfImages, 'double');
for ii = 1:numberOfImages
    idx = ceil(ii/imagesPerPosition);
    idx2 = randsample(size(barCarriers,3),1);
    images(:,:,ii) = barApertures(:,:,idx) .*  barCarriers(:,:,idx2);
end

% Recast as 8 bit integers
images = uint8(images + backgroundColor);


stimulus = [];
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.images     = images;
stimulus.seqtiming  = 0:1/frameRate:totalTime;
stimulus.seq        = zeros(size(stimulus.seqtiming))+blankImageIndex;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));


numIms =    round(stimulusDuration*frameRate);
numBlanks = round(isi*frameRate);
chunkSize = numIms + numBlanks;

for ii = 1:numberOfImages
    idx = (1:chunkSize) + (ii-1)*chunkSize;
    stimOrder = [repmat(ii, [1,numIms]) repmat(blankImageIndex, [1, numBlanks])];
    stimulus.seq(idx) = stimOrder;
end

% Add fixation sequence
minDurationInSeconds = 1;
maxDurationInSeconds = 5;
fixSeq = createFixationSequence(stimulus, 1/frameRate, minDurationInSeconds, maxDurationInSeconds);
stimulus.fixSeq = fixSeq;

maxUpdateInterval = 0.25;
stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);

% add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri' 
    otherwise
        stimulus.trigSeq  = double(stimulus.seq~=blankImageIndex);
end

% create stimulus.mat filename
fname = sprintf('%s_prf_%d.mat', stimParams.site, runNum);


% add table with elements to write to tsv file for BIDS
[~, Ai] = unique(stimulus.seq);

onset       = stimulus.seqtiming(Ai)';
duration    = ones(numberOfImages,1) * stimulusDuration;
trial_type  = ceil((1:numberOfImages)/imagesPerPosition)';
trial_name  = num2str(trial_type);
stim_file   = repmat(fname, numberOfImages,1);
stim_file_index = (1:numberOfImages)';

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% save
stimulus.display  = stimParams.display;
stimulus.modality = stimParams.modality;
stimulus.site     = stimParams.site;

save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

