function stimMakeFingermappingExperiment(stimParams, runNum,experimentType)

% Set a path to find .jpg files for now
resourcePath = fullfile(BAIRRootPath , 'motorStimuliResources', 'fingerMap');

%first load in the onsets file depending on the modality
if strcmpi(stimParams.modality, {'ecog', 'meg', 'eeg'})
    tmpData = load(fullfile(resourcePath,'stiminfo_ecog.txt'));
elseif strcmpi(stimParams.modality, 'fmri')
    tmpData = load(fullfile(resourcePath,'stiminfo_mri.txt'));
else
    error('Unknown Modality')
end

frameRate        = stimParams.display.frameRate;
onsets           = tmpData(:,1)/1000; %convert onsets from milliseconds to seconds
onsets           = round(onsets*frameRate)/frameRate; % match onset to frame rate
imgSeq           = tmpData(:,2);
experimentLength = max(onsets);

%initialize and set some stimulus properties
stimulus            = [];
stimulus.cat        = [20 21 22 23 24 25 26 27 28 29];
stimulus.categories = {'none', 'pinky', 'pinky:ring', 'pinky:middle','pinky:pointer',...
                        'all','thumb', 'thumb:pointer','thumb:middle','thumb:ring'};
stimulus.onsets     = onsets;
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.display    = stimParams.display;

stimulus.seqtiming  = 0:1/frameRate:experimentLength;
stimulus.fixSeq     = ones(size(stimulus.seqtiming));
stimulus.seq        = zeros(size(stimulus.seqtiming));

% Insert the image sequence into stimulus.seq
for ee = 1: length(onsets)-1
    startIdx = length(0:1/frameRate:(onsets(ee)));
    endIdx = length(0:1/frameRate:(onsets(ee+1)-1/frameRate));
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

% % Find the rectangle the image will be displayed in so we can shift the image into the center later
% screenRect       = size(zeros(stimulus.dstRect(4)-stimulus.dstRect(2): stimulus.dstRect(3)-stimulus.dstRect(1)));
% leftShift        = abs(.5*screenRect(2)-0.5*imgSize(2));
% topShift         = abs(.5*screenRect(1)-0.5*imgSize(1));
% shiftedLocation1 = topShift:topShift+imgSize(1)-1; %shifted to center
% shiftedLocation2 = leftShift:leftShift+imgSize(2)-1; %shifted to center

% Pre-allocate arrays to store images
images = zeros([size(tempImg) length(stimulus.categories)], 'uint8');

% % make a blank to insert between simulus presentations
% blankImg = images(:,:,:,1);
% blankImg(:) = 127;
 
% Load the images and resize them
for cc = 1:length(stimulus.cat)
    thisImgName = sprintf('%02d.png', cc);
    thisimage = imread(fullfile(imgFiles(cc).folder, thisImgName));
%     image = imresize(thisImage, [imgSize(1) imgSize(2)]);
%     images(:,:,:,cc) = 127; %first set the entire image to gray
    images(:,:,:,cc) = thisimage; %then insert the bitmap
    
end
% images(:,:,:,length(stimulus.cat)+1) = blankImg;
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
duration        = diff(onsets);
duration(end+1) = 1/frameRate; %the last onset is only on for a frame
trial_type      = stimulus.cat(imgSeq)';
trial_name      = stimulus.categories(imgSeq)';
stim_file       = repmat(fname, length(onset),1);
stim_file_index = repmat('n/a', length(stimulus.onsets),1);

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

% Save
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname))
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end
