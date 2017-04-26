function stimMakeTaskExperiment(s_example, modality)
% stimMakeTaskExperiment(s_example, modality)

switch(modality)
    case 'fMRI'
        numruns = 2;
    case 'MEG'
        numruns = 2;
end

for runnum = 1:numruns
    stimulus = [];
    stimulus.cmap       = s_example.stimulus.cmap;
    stimulus.srcRect    = s_example.stimulus.srcRect;
    stimulus.dstRect    = s_example.stimulus.dstRect;
    stimulus.images     = s_example.stimulus.images(:,:,end);
    stimulus.seqtiming  = (0:4*216-1)/4;
    stimulus.seq        = ones(size(stimulus.seqtiming));
    stimulus.seqtiming  = (0:4*216-1)/4;
    stimulus.seq        = ones(size(stimulus.seqtiming));
    stimulus.fixSeq     = ones(size(stimulus.seqtiming));
    
    idx = mod(stimulus.seqtiming,24)<12;
    stimulus.fixSeq(idx) = 3;
    
    % insert blips
    n = 54;
    mint = 1.5;
    maxt = 8;
    x = linspace(0,1,n) * (1-exp(-mint)) + exp(-mint);
    y= -log(x)/mint;
    y = y*(maxt-mint)+mint;
    
    y = round(y*4);
    stim_seq = randperm(n);
    y = y(stim_seq);
    blips = cumsum(y);
    
    stimulus.fixSeq(blips) = stimulus.fixSeq(blips)+1;
    
    switch lower(modality)
        case 'fmri'
        otherwise
            stimulus.diodeSeq = stimulus.fixSeq > 2;
            stimulus.trigSeq  = stimulus.fixSeq;
    end
    
    
    fname = sprintf('task_%s_%d', modality, runnum);
    save(fullfile(vistadispRootPath, 'Retinotopy', 'storedImagesMatrices', fname), 'stimulus')
    
end

