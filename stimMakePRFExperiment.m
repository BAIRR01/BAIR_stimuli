function stimMakePRFExperiment(s_example, modality)

imsize = s_example.stimulus.srcRect(4);
 
readPth  = 'https://wikis.nyu.edu/download/attachments/85394548/bar_apertures.mat?api=v2';
stimDir  = fullfile(BAIRRootPath, 'stimuli');
fname    = 'bar_apertures.mat';
writePth = fullfile(stimDir, fname);
websave(writePth,readPth);
im = load(writePth);

bar_apertures = imresize(im.bar_apertures, imsize, 'nearest');
bar_carrier   = stimMakeBarCarrier();


numim = size(bar_apertures,3)*3;

% 3 images per aperture
images = zeros(imsize, imsize, numim, 'uint8')+127;
for ii = 1:numim
    idx = ceil(ii/3);
    idx2 = randsample(size(bar_carrier,3),1);
    images(:,:,ii) = bar_apertures(:,:,idx) .* (bar_carrier(:,:,idx2)-.5) * 255+127;
end

switch(modality)
    case 'fMRI'
        numruns = 4;
    case 'MEG'
        numruns = 12;
end
for runnum = 1:numruns
    stimulus = [];
    stimulus.cmap       = s_example.stimulus.cmap;
    stimulus.srcRect    = s_example.stimulus.srcRect;
    stimulus.dstRect    = s_example.stimulus.dstRect;
    stimulus.images     = images;
    stimulus.seqtiming  = (0:numim-1)/3;
    stimulus.seq        = 1:length(stimulus.seqtiming);
    stimulus.fixSeq     = ones(size(stimulus.seqtiming));
    
    switch lower(modality)
        case 'fmri'
        otherwise
            stimulus.trigSeq  = double(stimulus.seq>0);
            stimulus.diodeSeq = stimulus.trigSeq;
    end
    
    
    fname = sprintf('ret_%s_%d', modality, runnum);
    save(fullfile(vistadispRootPath, 'Retinotopy', 'storedImagesMatrices', fname), 'stimulus')
    
end
