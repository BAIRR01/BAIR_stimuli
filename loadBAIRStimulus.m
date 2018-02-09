function [master_stimulus] = loadBAIRStimulus(stimulusType, site, runNum)

% LOAD MASTER STIMULUS 

readPth = vistadispRootPath; 

% Source needs to be updated to location where stimuli are stored and shared, e.g. Flywheel
% use SciTran software: https://github.com/scitran/client, create toolboxes
% to add to bairStimuli json files
%
% older example code for downloading files:
% readPth = sprintf('https://wikis.nyu.edu/download/attachments/85394548/bar_carrier%d.mat?api=v2', ii);
% stimDir = fullfile(BAIRRootPath, 'stimuli');
% fname = sprintf('bar_carrier%d.mat', ii);
% writePth = fullfile(stimDir, fname);
% if ~exist(writePth, 'file'), websave(writePth,readPth); end
% im = load(writePth);

fname = sprintf('%s_%s_%d', stimulusType, site, runNum);
load(fullfile(readPth, 'Retinotopy', 'storedImagesMatrices',  fname), 'stimulus')

master_stimulus= stimulus;

end
