function stimMakePRFExperiment(stimParams, runNum, stimulusDuration,dwellTimePerImage, isi)

% Load aperture images
barApertures = stimMakePRFApertures(stimParams);

% Load carrier images
barCarriers   = stimMakePRFCarriers(stimParams);

% How many images shown for each aperture?
imagesPerPosition = round(1 / (stimulusDuration + isi));

% Total number of images
numberOfImages = imagesPerPosition * size(barApertures,3);

% First image is a blank
blankImageIndex  = 1;

% Total time points
numberOfTimePoints = 1/dwellTimePerImage * size(barApertures,3);

% Image size
imageSizeInPixels = size(stimParams.stimulus.images);

% Recast to double centered at background color
barCarriers = double(barCarriers);
backgroundColor = mode(barCarriers(:));
barCarriers = barCarriers  - backgroundColor;

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
stimulus.seqtiming  = (0:numberOfTimePoints-1)*dwellTimePerImage;
stimulus.seq        = zeros(size(stimulus.seqtiming))+blankImageIndex;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));


numIms =    round(stimulusDuration/dwellTimePerImage);
numBlanks = round(isi/dwellTimePerImage);
chunkSize = numIms + numBlanks;

for ii = 1:numberOfImages
    idx = (1:chunkSize) + (ii-1)*chunkSize;
    stimOrder = [repmat(ii, [1,numIms]) repmat(blankImageIndex, [1, numBlanks])];
    stimulus.seq(idx) = stimOrder;
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
        stimulus.trigSeq  = double(stimulus.seq~=blankImageIndex);
        stimulus.diodeSeq = stimulus.trigSeq;
end

% create stimulus.mat filename
fname = sprintf('prf_%s_%d', stimParams.site, runNum);

% add table with elements to write to tsv file for BIDS
onset       = stimulus.seqtiming((0:numberOfImages-1)*chunkSize+1)';
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

