function [masterImages] = loadBAIRStimulus(stimulusType, site, runNum)

% LOAD MASTER STIMULUS 

readPth = vistadispRootPath; 

% Source needs to be updated to location where stimuli are stored and shared, e.g. Flywheel
% example code:
% readPth = sprintf('https://wikis.nyu.edu/download/attachments/85394548/bar_carrier%d.mat?api=v2', ii);
% stimDir = fullfile(BAIRRootPath, 'stimuli');
% fname = sprintf('bar_carrier%d.mat', ii);
% writePth = fullfile(stimDir, fname);
% if ~exist(writePth, 'file'), websave(writePth,readPth); end
% im = load(writePth);

fname = sprintf('hrf%s_%s_%d', stimulusType, site, runNum);
load(fullfile(readPth, 'Retinotopy', 'storedImagesMatrices',  fname), 'stimulus')

masterImages = stimulus.images;

end
