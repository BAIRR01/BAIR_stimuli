function barApertures = stimMakePRFApertures(stimParams)

imageSizeInPixels = size(stimParams.stimulus.images);

stimDir  = fullfile(BAIRRootPath, 'stimuli');
fname    = 'bar_apertures.mat';
writePth = fullfile(stimDir, fname);
if ~exist(writePth, 'file')
    readPth  = 'https://wikis.nyu.edu/download/attachments/85394548/bar_apertures.mat?api=v2';
    websave(writePth,readPth);
end
im = load(writePth);

barApertures = imresize(im.bar_apertures, imageSizeInPixels, 'nearest');

return