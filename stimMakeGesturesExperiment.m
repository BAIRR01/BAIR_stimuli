function stimMakeGesturesExperiment(stimParams,  runNum, TR, stimDurationSeconds,experimentType)
%
% NOTE: This is just a framework, will update soon
% Commented text is in need of updating or checking
%

% Set a path to find .txt and .jpg files for now
resourcePath = fullfile(BAIRRootPath , 'motorStimuliResources');

% contains for each bitmap, the fMRI pulse count on which it should be shown
rawOnsets             = load('picture_onset_sequence.txt');
onsets                = round(rawOnsets/TR)*TR;
numberofEvents        = length(onsets);
preScanPeriod         = round(12/TR)*TR;
postScanPeriod        = preScanPeriod;
eventLength           = round(stimDurationSeconds/TR)*TR;
experimentLength      = max(onsets) + eventLength + preScanPeriod + postScanPeriod;

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

stimulus = [];
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;
stimulus.onsets     = onsets;
stimulus.cat        = [1 2 3 4];
stimulus.categories = { 'D', 'F', 'V', 'Y'};

stimulus.seqtiming  = 0:1/frameRate:experimentLength;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));
stimulus.seq        = zeros(size(stimulus.seqtiming)); %initialize it for now

% contains the filenames of the bitmaps to be shown
testImgSeq     = importdata('bitmap_filename_sequence.txt');
imgLetterSeq   = string(zeros(length(testImgSeq),1));
imgNumSeq      = zeros(length(testImgSeq),1);
for ii = 1:length(stimulus.cat)
    idx               = contains(testImgSeq, sprintf('exec_stim_%d', ii));
    imgNumSeq(idx)    = stimulus.cat(ii);
    imgLetterSeq(idx) = stimulus.categories{ii};
end

eventLengthInFrames = length(0:1/frameRate:eventLength);
for ee = 1: numberofEvents
    StartFrame = length(0:1/frameRate:preScanPeriod+onsets(ee));
    EndFrame = StartFrame + eventLengthInFrames;
    stimulus.seq(StartFrame:EndFrame) = imgNumSeq(ee);
end

% first, find all the bitmaps
bitmapPth = fullfile(resourcePath, 'bitmaps');
files     = dir([bitmapPth '/*jpg']);

switch experimentType
    case 'GESTURESTRAINING'
        %figure out which ones are for training
        trainingIdx = contains({files.name},stimulus.categories);
        imgFiles    = files(trainingIdx);
    case 'GESTURES'
        %figure out which ones are for testing
        testIdx  = contains({files.name},'exec_stim');
        imgFiles = files(testIdx);
end

%load in the first bitmap to set a standard size for the rest in case they
% are not the same (this is the case with the training images)

%imageSizeInPixels = size(stimParams.stimulus.images); % based on visual simuli
[tempImg, ~,~] = imread(fullfile(imgFiles(1).folder, imgFiles(1).name));
%resizedImgSize = size(imresize(tempImg, [imageSizeInPixels(1) NaN]));
imgSize = size(tempImg);

% Pre-allocate arrays to store images
images = zeros([imgSize length(stimulus.categories)], 'uint8');

% Load the images and resize them
for cc = 1:length(stimulus.cat)
   %check and see that this order matches
    imageForThisTrial = imread(fullfile(files(cc).folder, files(cc).name));
    image = imresize(imageForThisTrial, [imgSize(1) imgSize(2)]);
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

% Create stim_file name
fname = sprintf('%s_%s_%d.mat', stimParams.site,lower(experimentType), runNum);

% Add table with elements to write to tsv file for BIDS
onset           = onsets;
duration        = ones(size(onset)) * eventLength;
trial_type      = imgLetterSeq;
trial_name      = imgLetterSeq;
stim_file       = repmat(fname, length(onset),1);
stim_file_index = repmat('n/a', length(stimulus.onsets),1);

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% Save
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname))
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end

