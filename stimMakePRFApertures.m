function barApertures = stimMakePRFApertures(stimParams)

imageSizeInPixels = size(stimParams.stimulus.images);

% readPth  = 'https://wikis.nyu.edu/download/attachments/85394548/bar_apertures.mat?api=v2';
% stimDir  = fullfile(BAIRRootPath, 'stimuli');
% fname    = 'bar_apertures.mat';
% writePth = fullfile(stimDir, fname);
% websave(writePth,readPth);
% im = load(writePth);

load('~/Desktop/bar_apertures.mat'); 
barApertures = imresize(im.bar_apertures, imageSizeInPixels, 'nearest');

return