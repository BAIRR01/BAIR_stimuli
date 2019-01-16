function stimMakeGesturesExperiment(stimParams,  runNum, TR, stimDurationSeconds,experimentType)
% stimMakeGesturesExperiment(stimParams,  runNum, TR, eventLength,experimentType)
%
% Loads in bitmaps found in ~/BAIRstimuli/motorStimuliResources/bitmaps
%
% This code should be run from s_Make_BAIR_Motor_Experiments

% Set a path to find .jpg files for now
resourcePath = fullfile(BAIRRootPath , 'motorStimuliResources');

stimulus = [];
stimulus.cat        = [10 11 12 13];
stimulus.categories = { 'D', 'F', 'V', 'Y'};

switch(lower(stimParams.modality))
    case 'fmri'
        if strcmp(experimentType,'GESTURES')
            isiMin         = round(6/TR)*TR; % 6 seconds
            isiMax         = round(15/TR)*TR; % 15 seconds
            preScanPeriod  = round(12/TR)*TR; % 12 seconds
            desiredLength  = round(480/TR)*TR; % 480 seconds (8 min)
            stimDurationSeconds    = round(stimDurationSeconds/TR)*TR;
            
        elseif contains(experimentType,{'GESTURESPRACTICE','GESTURESLEARNING'})
            % keep these shorter and with the same timing as ecog/meg
            isiMin        = 4; % seconds
            isiMax        = 6; % seconds
            preScanPeriod = 3; % seconds
            desiredLength = 300; % seconds
        end
        
    case {'ecog' 'eeg' 'meg'}
        isiMin        = 4; % seconds
        isiMax        = 6; % seconds
        preScanPeriod = 3; % seconds
        desiredLength = 480; % seconds
    otherwise
        error('Unknown modality')
end

rng('shuffle');
% Find the number of events we can fit in desired experiment time and jitter ISIs
numberofEvents = floor(desiredLength / (stimDurationSeconds+((isiMax + isiMin) / 2)));
possibleISIs   = linspace(isiMin,isiMax,numberofEvents-1);
isiSeq         = randperm(numberofEvents-1);

% Find the onsets and match the stimulus presentation to the frame rate
frameRate        = stimParams.display.frameRate;
onsets           = cumsum([preScanPeriod stimDurationSeconds+possibleISIs(isiSeq)]);
onsets           = round(onsets*frameRate)/frameRate;
onsetFrameIdx    = round(onsets*(frameRate*.5));
postScanPeriod   = preScanPeriod;
experimentLength = max(onsets) + stimDurationSeconds + postScanPeriod;

% set some stimulus properties
stimulus.ISI        = possibleISIs;
stimulus.prescan    = preScanPeriod;
stimulus.postscan   = postScanPeriod;
stimulus.onsets     = onsets;
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;

stimulus.seqtiming  = 0:(1/frameRate)*2:experimentLength;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));
stimulus.seq        = zeros(size(stimulus.seqtiming));

% Figure out a random order to present the images
imgSeq      = randi([1,length(stimulus.cat)],length(onsets),1);

eventLengthInFrames = round(length(0:(1/frameRate)*2:stimDurationSeconds));
for ee = 1: numberofEvents
    stimulus.seq (onsetFrameIdx(ee):onsetFrameIdx(ee)+eventLengthInFrames) = imgSeq(ee);
    stimulus.fixSeq(onsetFrameIdx(ee):onsetFrameIdx(ee)+eventLengthInFrames) = 0;
end
blankIdx = stimulus.seq == 0;
stimulus.seq(blankIdx) = length(stimulus.cat)+1;


switch experimentType
    case 'GESTURESLEARNING'
        % first, find all the bitmaps
        bitmapPth = fullfile(resourcePath,'gestures', 'bitmaps');
        files     = dir([bitmapPth '/*jpg']);
        
        %figure out which ones are for training
        trainingIdx = contains({files.name},stimulus.categories);
        imgFiles    = files(trainingIdx);
        
        %load in one image to resize and use that for the experiment
        [tempImg, ~,~] = imread(fullfile(imgFiles(1).folder, imgFiles(1).name));
        imgSize = size(imresize(tempImg, 2));
        isLoadedImg = 1;
        
    case {'GESTURESPRACTICE' ,'GESTURES'}
        % load in letters from Kay et al, 2013 PLOS CB and find the ones we need
        load (fullfile(resourcePath,'gestures', 'letters.mat'), 'letters');
        
        %the corresponding letter numbers for our stimuli
        lettersIdx = [4 6 22 25];
        gestureLetters = letters(:,:,lettersIdx);
        
        % set the colors equal to what we're using for this experiment
        tmpIdx = gestureLetters == 0;
        gestureLetters (tmpIdx) = 127; % make the background gray
        tmpIdx = gestureLetters < 0;
        gestureLetters (tmpIdx) = 0; % make the letter black
        
        % there is a lot of space around the letters, so crop it
        origImgSize = size(gestureLetters);
        cropAmt = round([(origImgSize(1))/2 (origImgSize(1))/2]);
        cropIdx1 = cropAmt(1)/2:origImgSize(1)-cropAmt(1)/2;
        cropIdx2 = cropAmt(2)/2:origImgSize(2)-cropAmt(2)/2;
        gestureLetters = gestureLetters(cropIdx1,cropIdx2,:);
        
        % we want to make the letters a little bigger, so get the size
        imgSize = size(imresize(gestureLetters(:,:,1), 3));
        isLoadedImg = 0;
end

% Find the rectangle the image will be displayed in so we can shift the image into the center later
screenRect       = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));
leftShift        = round(abs(.5*screenRect(2)-0.5*imgSize(2)));
topShift         = round(abs(.5*screenRect(1)-0.5*imgSize(1)));
shiftedLocation1 = topShift:topShift+imgSize(1)-1; %shifted to center
shiftedLocation2 = leftShift:leftShift+imgSize(2)-1; %shifted to center

% Pre-allocate arrays to store images
images = zeros([screenRect 3 length(stimulus.categories)+1], 'uint8');

% make a blank to insert between simulus presentations
blankImg = images(:,:,:,1);
blankImg(:) = 127;

% Load the images and resize them
for cc = 1:length(stimulus.cat)
    %first set the entire image to gray
    images(:,:,:,cc) = 127;
    
    if isLoadedImg
        thisImage = imread(fullfile(imgFiles(cc).folder, imgFiles(cc).name));
    else
        thisImage = ones([size(gestureLetters(:,:,cc)) 3]);
        for jj = 1:3
            thisImage(:,:,jj) = gestureLetters(:,:,cc);
        end
    end
    % then insert the resized image in the center
    images(shiftedLocation1,shiftedLocation2,:,cc) = imresize(thisImage, imgSize(:,:));
end
images(:,:,:,length(stimulus.cat)+1) = blankImg;
stimulus.images     = images;

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % Write trigger sequence
        stimulus.trigSeq        = zeros(length(stimulus.seq),1);
        idx                     = find(diff(stimulus.seq < max(stimulus.seq)) == -1);
        stimulus.trigSeq(idx)   =  stimulus.cat(stimulus.seq(idx));
        stimulus.trigSeq(1)     = 255; %experiment onset
        stimulus.trigSeq(end)   = 255; %experiment offset
end

% Create stim_file name
fname = sprintf('%s_%s_%d.mat', stimParams.site,lower(experimentType), runNum);

% Add table with elements to write to tsv file for BIDS
onset           = round(stimulus.onsets,3)';
duration        = repmat(round(stimDurationSeconds,3),length(onsets),1);
trial_type      = stimulus.cat(imgSeq)';
trial_name      = stimulus.categories(imgSeq)';
stim_file       = repmat(fname, length(onset),1);
stim_file_index = repmat('n/a', length(stimulus.onsets),1);

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% Save
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname))
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end

