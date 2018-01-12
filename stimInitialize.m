function stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg)
% Get the data structure for an example stimulus. The fields will be
% populated for actual stimuli.
%
% stimParams = stimInitialize(experimentSpecs, 1, 16.6);

displayParameters       = loadDisplayParams(experimentSpecs.displays{whichSite});

screenHeightInPixels    = displayParameters.numPixels(2);

screenWidthInPixels     = displayParameters.numPixels(1);

screenHeightInDegrees   = pix2angle(displayParameters, screenHeightInPixels);

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


stimParams.stimulus = stimulus;
stimParams.display  = displayParameters;
stimParams.modality = experimentSpecs.modalities(whichSite);
stimParams.experimentSpecs = experimentSpecs(whichSite,:);

return
