function stimMakeTaskExperiment(stimParams,  runNum, TR)
% stimMakeTaskExperiment(stimParams,  runNum, TR)

site = stimParams.experimentSpecs.sites{1};
% 9 alternations of 24 s (12 s task, 12 s no task)

numBlocks           = 18;
desiredBblockLength = 12; % seconds
trsPerBlock         = round(desiredBblockLength / TR);
blockLength         = trsPerBlock * TR;
experimentLength    = blockLength * numBlocks;
frameRate           = stimParams.display.frameRate;

stimulus = [];
stimulus.cmap       = stimParams.stimulus.cmap;
stimulus.srcRect    = stimParams.stimulus.srcRect;
stimulus.dstRect    = stimParams.stimulus.destRect;
stimulus.images     = stimParams.stimulus.images(:,:,end);
stimulus.seqtiming  = 0:1/frameRate:experimentLength;
stimulus.seq        = ones(size(stimulus.seqtiming));
stimulus.fixSeq     = ones(size(stimulus.seqtiming));





idx = mod(stimulus.seqtiming,blockLength*2)<blockLength;
stimulus.fixSeq(idx) = 3;

% insert blips
minISI    = 1.5; % seconds
maxISI    = 8;


fixSeq = createFixationSequence(stimulus, 1/frameRate, minISI, maxISI);
blips  = abs(diff([1 fixSeq]))>0;

stimulus.fixSeq(blips) = stimulus.fixSeq(blips)+1;


% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri' 
        % no trigger sequence needed
    otherwise
        % Write binary trigger sequence:
        stimulus.trigSeq   = blips;
end


maxUpdateInterval = 0.25;
stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);


% Create stim_file name
fname = sprintf('%s_dottask_%d.mat', site, runNum);


% Add table with elements to write to tsv file for BIDS
conditions  = {'task' 'rest'};
    
onset       = ((0:numBlocks-1)*blockLength)';
duration    = ones(size(onset)) * blockLength;
trial_type  = (mod(0:numBlocks-1,2)+1)';
trial_name  = conditions(trial_type)';
stim_file   = repmat(fname, length(onset),1);
stim_file_index = ones(size(onset));

stimulus.tsv = table(onset, duration, trial_type, trial_name, stim_file, stim_file_index);

save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus')

end


