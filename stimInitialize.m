function s_example = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg)
% Get the data structure for an example stimulus. The fields will be
% populated for actual stimuli.
%
% s_example = stimInitialize(experimentSpecs, 1, 16.6);

displayParameters       = loadDisplayParams(experimentSpecs.displays{whichSite});

screenHeightInPixels    = displayParameters.numPixels(2);

screenWidthInPixels     = displayParameters.numPixels(1);

screenHeightInDegrees   = 2*pix2angle(displayParameters, screenHeightInPixels/2);

fractionOfScreenToUse   = stimDiameterDeg/screenHeightInDegrees;

stimSizeInPixels        = fractionOfScreenToUse * screenHeightInPixels;

% We want our image size to be an even number of pixels for subsequent
% image processing
nearestEven = @(x) round(x/2)*2;
stimSizeInPixels = nearestEven(stimSizeInPixels);



% colormap: 1:256 (gamma table is controlled separately)
stimulus.cmap    = (ones(3,1)*(1:256))'; 

% image size and position
stimulus.srcRect = [0 0 stimSizeInPixels stimSizeInPixels];
left   = (screenWidthInPixels - stimSizeInPixels)/2;
top    = (screenHeightInPixels-stimSizeInPixels)/2;
right  = screenWidthInPixels - left; 
bottom = screenHeightInPixels - top;

stimulus.destRect = [left, top, right, bottom];

% a blank image
stimulus.images = uint8(zeros(stimSizeInPixels)+127);

stimulus.seqtiming  = [];
stimulus.seq        = [];
stimulus.fixSeq     = [];


s_example.stimulus = stimulus;
s_example.display  = displayParameters;
s_example.modality = experimentSpecs.modalities(whichSite);

return

% % load example file, stored on the web
% readPth  = 'https://wikis.nyu.edu/download/attachments/85394548/BAIRstimExample.mat?api=v2';
% 
% % local directory where it will be copied
% stimDir  = fullfile(BAIRRootPath, 'stimuli');
% fname    = 'BAIRstimExample.mat';
% writePth = fullfile(stimDir, fname);
% 
% % check whether local directory exists
% if ~exist(stimDir, 'dir'), mkdir(stimDir); end
% addpath(stimDir);
% 
% % copy from web to local
% websave(writePth,readPth);
% 
% % load it into the workspace
% s_example = load(writePth);
