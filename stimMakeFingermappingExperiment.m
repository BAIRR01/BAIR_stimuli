function stimMakeFingermappingExperiment(stimParams, runNum,experimentType)

% Set a path to find .jpg files for now
resourcePath = fullfile(BAIRRootPath , 'motorStimuliResources', 'fingerMap');

%first load in the onsets file depending on the modality
if any(strcmpi(stimParams.modality, {'ecog', 'meg', 'eeg'}))
    tmpData = load(fullfile(resourcePath,'stiminfo_ecog.txt'));
elseif strcmpi(stimParams.modality, 'fmri')
    tmpData = load(fullfile(resourcePath,'stiminfo_mri.txt'));
else
    error('Unknown Modality')
end

frameRate        = stimParams.display.frameRate;
onsets           = tmpData(:,1)'/1000; %convert onsets from milliseconds to seconds
onsets           = round(onsets*frameRate)/frameRate; % match onset to frame rate
imgSeq           = tmpData(:,2);
experimentLength = max(onsets);

%initialize and set some stimulus properties
stimulus            = [];
stimulus.cat        = [20 21 22 23 24 25 26 27 28 29];
stimulus.categories = {'open palm', 'pinky close', 'ring close', 'middle close','pointer close',...
                        'closed palm','thumb open', 'pointer open','middle open','ring open'};
stimulus.onsets     = onsets;
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;

stimulus.seqtiming  = 0:((1/frameRate)*2):experimentLength;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));
stimulus.seq        = zeros(size(stimulus.seqtiming));

% Insert the image sequence into stimulus.seq
for ee = 1: length(onsets)-1
    startIdx = length(0:(1/frameRate)*2:(onsets(ee)));
    endIdx = length(0:(1/frameRate)*2:(onsets(ee+1)-(1/frameRate)*2));
     stimulus.seq(startIdx:endIdx) = imgSeq(ee);
end
% Insert the last stimulus.seq value we missed
stimulus.seq(end) = imgSeq(end);

% Figure out which images to use based on which hand we are using
switch experimentType
    case 'FINGERMAPPINGLEFT'
        imgDir = fullfile(resourcePath, '02_ExeSeq_pics', '*.png');
    case 'FINGERMAPPINGRIGHT'
        imgDir = fullfile(resourcePath, '01_ExeSeq_pics', '*.png');
end
% Find the images
imgFiles = dir(imgDir);

%load in one image to set size for the experiment
[tempImg, ~,~] = imread(fullfile(imgFiles(1).folder, imgFiles(1).name));
imgSize = size(tempImg);

% Find the destination rectangle so we can crop our stimuli
screenRect  = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));

% Pre-allocate arrays to store images
images = zeros([screenRect imgSize(3) length(stimulus.categories)], 'uint8');
% we're going to gray out the center fixation point, so find the image center
imgCenter = 0.5*imgSize;

cropAmt    = 500; 
cropIdx1 = round(imgCenter(1)-cropAmt:imgCenter(1)+cropAmt);
cropIdx2 = round(imgCenter(2)-cropAmt:imgCenter(2)+cropAmt);

% Load the images and resize them
for cc = 1:length(stimulus.cat)
    thisImgName = sprintf('%02d.png', cc);
    thisImg = imread(fullfile(imgFiles(cc).folder, thisImgName));
    %the gray in these images is off, set them equal to a different gray
    grayIdx = thisImg == 148;  
    thisImg(grayIdx) = 127;
    %the fixation point is in the center, so make everything in the center gray 
    thisImg(imgCenter(1)-10:imgCenter(1)+10,imgCenter(2)-10:imgCenter(2)+10, :) = 127;
    images(:,:,:,cc) =  imresize(thisImg(cropIdx1,cropIdx2, :), screenRect);
    
end
stimulus.images     = images;

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri'
        % no trigger sequence needed
    otherwise
        % Write trigger sequence
        stimulus.trigSeq        = zeros(length(stimulus.seq),1);
        idx                     = find(diff(stimulus.seq) ~= 0);
        stimulus.trigSeq(idx)   =  stimulus.cat(stimulus.seq((idx+1)));
        stimulus.trigSeq(1)     = 255; %experiment onset
        stimulus.trigSeq(2)     = stimulus.cat(stimulus.seq((idx(1)))); %first image
        stimulus.trigSeq(end)   = 255; %experiment offset
end

% Create stim_file name
fname = sprintf('%s_%s_%d.mat', stimParams.site,lower(experimentType), runNum);

% Add table with elements to write to tsv file for BIDS
onset           = round(stimulus.onsets(1:end-1),3);
duration        = round(diff(onsets),3);
trial_type      = stimulus.cat(imgSeq(1:end-1))';
trial_name      = stimulus.categories(imgSeq(1:end-1))';
stim_file       = repmat(fname, length(onset),1);
stim_file_index = repmat('n/a', length(onset),1);

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% Save
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname))
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end
