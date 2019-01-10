function stimMakeGesturesExperiment(stimParams,  runNum, TR, stimDurationSeconds,experimentType)
%
% NOTE: This is just a framework, will update soon
% Commented text is in need of updating or checking
%

% Set a path to find .txt and .jpg files for now
resourcePath = fullfile(BAIRRootPath , 'motorStimuliResources');

% contains for each bitmap, the fMRI pulse count on which it should be shown
onsets = load('picture_onset_sequence.txt');
% contains the filenames of the bitmaps to be shown
imgSeq = importdata('bitmap_filename_sequence.txt');

% Initialize, then set some the stimulus parameters
frameRate           = stimParams.display.frameRate;

% from provided UMCU code (generate_stimuli.m) 
% scan_period = 850; % for mri - equal to Presentation equivalent 
% isi_min =  6; % minimum inter-stimulus distance, between onset times, in s
% isi_max = 15; % maximum inter-stimulus distance, between onset times, in s
% max_dur = 480; % maximum duration of the whole task (seconds)
% isi_min1 = round(isi_min / scan_period* 1000);
% isi_max1 = round(isi_max / scan_period* 1000);
% max_dur1 = round(max_dur / scan_period* 1000);
% n_events = floor(max_dur / ((isi_max + isi_min) / 2));
% isi = [1; randi(isi_max - isi_min + 1, n_events - 1, 1) + isi_min - 1];
% onsets = cumsum(isi);
% events = randi(length(stimulus.cat), n_events, 1);

% this is from stimMakeTaskExperiment and needs to be updated
% numBlocks           = 18;
% desiredBblockLength = 12; % seconds
% trsPerBlock         = round(desiredBblockLength / TR);
% blockLength         = trsPerBlock * TR;
% experimentLength    = blockLength * numBlocks;

stimulus = [];
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;

stimulus.cat        = [1 2 3 4];
stimulus.categories = { 'D', 'F', 'V', 'Y'};

% This needs to get fixed
% stimulus.seqtiming  = 0:1/frameRate:experimentLength;
% stimulus.seq        = ones(size(stimulus.seqtiming));
% stimulus.fixSeq     = ones(size(stimulus.seqtiming));

% first, find all the bitmaps
bitmapPth = fullfile(resourcePath, 'bitmaps');
files = dir([bitmapPth '/*jpg']);

switch experimentType
    case 'GESTURESTRAINING'
        %figure out which ones are for training
        trainingIdx = contains({files.name},stimulus.categories);
        imgFiles = files(trainingIdx);
    case 'GESTURES'
        %figure out which ones are for testing
        testIdx = contains({files.name},'exec_stim');
        imgFiles = files(testingIdx);
end

%load in the first bitmap to set a standard size for the rest in case they
% are not the same (this is the case with the training images)
imageSizeInPixels = size(stimParams.stimulus.images); % based on visual simuli
tempImg = imread(fullfile(imgFiles(1).folder, imgFiles(1).name));
%resize
resizedImgSize = size(imresize(tempImg, [imageSizeInPixels(1) NaN]));

% Pre-allocate arrays to store images
images = zeros([resizedImgSize length(stimulus.categories)], 'uint8');

% Create the stimuli
for cc = 1:length(stimulus.cat)
    whichImg =  imgFiles(cc);
    imageForThisTrial = imread(fullfile(files(whichImg).folder, files(whichImg).name));
    image = imresize(imageForThisTrial, [resizedImgSize(1) resizedImgSize(2)]);
    images(:,:,:,cc) = image;
end

stimulus.images     = images;

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % Write binary trigger sequence:
        stimulus.trigSeq   = blips;
end

% % Sparsify stimulus sequence
% maxUpdateInterval = 0.25;
% stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);

% Create stim_file name
fname = sprintf('%s_%s_%d.mat', site,lower(experimentType), runNum);

% Add table with elements to write to tsv file for BIDS
conditions  = { 'D', 'F', 'V', 'Y'};
% condition_idx = mod(0:numBlocks-1,2)+1;

% onset       = ((0:numBlocks-1)*blockLength)';
% duration    = ones(size(onset)) * blockLength;
% trial_type  = stimulus.cat(condition_idx)';
% trial_name  = conditions(condition_idx)';
% stim_file   = repmat(fname, length(onset),1);
% stim_file_index = ones(size(onset));

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% Save
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end

