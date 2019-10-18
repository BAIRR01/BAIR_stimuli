function stimMakeLocalizerExperiment(stimParams, runNumber, stimulusType, onsetTimeMultiple, TR)

%% LOCALIZER EXPERIMENT

%% Make the images

% Determine if we're creating the master or loading & resizing for a specific display
site = stimParams.experimentSpecs.sites{1};
imageSizeInPixels = size(stimParams.stimulus.images);

switch stimulusType
    case 'SIXCATLOC'

        fprintf('[%s]: Creating stimulus file for stimulusType: SIXCATLOC, runID: %d.\n',mfilename, runNumber);

        categories = {...
            'bodies' ...
            'buildings' ...
            'faces' ...
            'objects' ...
            'scenes' ...
            'scrambled' ...
            };

        categoryNumberToAdd  = 0; % to make sure we have UNIQUE category number across all localizers
        numberOfImagesPerCat = 24;

        % Pre-allocate arrays to store images
        images = zeros([imageSizeInPixels  3 length(categories) * numberOfImagesPerCat], 'uint8');
        im_cell = cell([1 length(categories)]);
        catindex = zeros(1, length(categories) * numberOfImagesPerCat);
        imCount = 1;

        % Download original, unfiltered face, house, letter stimuli
        fprintf('[%s]: Loading stimuli...\n',mfilename);

        stimDir  = fullfile(BAIRRootPath, 'stimuli');
        fname    = 'sixcatlocalizer.mat';
        load(fullfile(stimDir, fname));

        % Category-specific settings
        numberOfCategories = length(categories);
        
                
        % Create the stimuli
        for cc = 1:numberOfCategories

            fprintf('[%s]: Creating stimuli at %d x %d pixels resolution: %s.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2), categories{cc});
            imageArray = eval(categories{cc});    
            totalNumberOfImagesAvailable = size(imageArray,4);
            % Remove a few stimuli with bad backgrounds
            switch categories{cc}
                case 'objects'
                    exclude = [121 131 141 159 172]; % These are natural objects
                    exclude = [exclude exclude-120]; % also remove same number of manmade objects
                    ind = setdiff(1:totalNumberOfImagesAvailable, exclude);
                case 'bodies'
                    exclude1 = [1 12]; % feet
                    exclude2 = [126 131 132 134 140 144 146 148 152:155 161 165 168 172 173 204 207 210 214 221 225 240]; % hands
                    exclude = [exclude1 exclude2(1:(end-length(exclude1)))-120]; % also remove same number feet
                    ind = setdiff(1:totalNumberOfImagesAvailable, exclude);
                otherwise
                    ind = 1:totalNumberOfImagesAvailable;
            end
            imageArray = imageArray(:,:,:,ind);            
            totalNumberOfImagesAvailable = length(ind);
                    
            % Pick which stimuli to select from original set
            % ODD for runNum == 1, EVEN for runNum == 2;
            if mod(runNumber,2) ~= 0
                startInx = 1;
            else
                startInx = 2;
            end
            switch categories{cc}
                case {'bodies', 'faces', 'objects'}
                    % bodies: 1:120 are feet, 121:240 are hands
                    % faces: 1:120 are female, 121:240 are male
                    % objects: 1:120 are manmade, 121:240 are natural
                    numberOfImagesPerSubCat = numberOfImagesPerCat/2;
                    totalNumberOfImagesAvailable = totalNumberOfImagesAvailable/2;
                    imageIndex = startInx:totalNumberOfImagesAvailable/numberOfImagesPerSubCat:totalNumberOfImagesAvailable;
                    imageIndex = [imageIndex imageIndex+totalNumberOfImagesAvailable];
                case {'buildings', 'scrambled'}
                    imageIndex = startInx:totalNumberOfImagesAvailable/numberOfImagesPerCat:totalNumberOfImagesAvailable;
                case 'scenes'
                    % scenes: 1:80 are indoor, 81:160 are outdoor
                    % manmade, 161:240 are outdoor natural
                    numberOfImagesPerSubCat = round(numberOfImagesPerCat/3);
                    totalNumberOfImagesAvailable = totalNumberOfImagesAvailable/3;
                    imageIndex1 = startInx:totalNumberOfImagesAvailable/numberOfImagesPerSubCat:totalNumberOfImagesAvailable;
                    imageIndex = [imageIndex1 imageIndex1+totalNumberOfImagesAvailable imageIndex1+totalNumberOfImagesAvailable*2];
                    imageIndex = round(imageIndex);
            end
            imageIndex = round(imageIndex);

            for ii = 1:numberOfImagesPerCat
                %inputImage = imageArray(:,:,imageIndex(ii));
                %imageForThisTrial = createFilteredStimulus(stimParams,inputImage,imageProcessingParams);

                % Double to unsigned 8 bit integer, needed for vistadisp
                %image8Bit = uint8((imageForThisTrial+.5)*255);
                %images(:,:,imCount) = image8Bit;
                %im_cell{categoryIndex(cc)}(:,:,ii) = image8Bit;

                inputImage = imageArray(:,:,:,imageIndex(ii));
                inputImage = imresize(inputImage, imageSizeInPixels);

                images(:,:,:,imCount) = inputImage;
                im_cell{cc}(:,:,:,ii) = inputImage;
                catindex(imCount) = cc+categoryNumberToAdd;
                imCount = imCount + 1;
            end
        end
        
        % Make sure images that contain grayscale pixels match the background
        backgroundColor = mode(images(:));
        fprintf('[%s]: Fixing stimulus backgrounds...\n',mfilename);
        for ii = 1:size(images,4)
            if ~contains(categories(catindex(ii)), 'scenes') % don't do this for the scenes
                I = images(:,:,:,ii);
                Imode = mode(I(:));
                if Imode ~= backgroundColor
                   ind = sum(I==Imode,3)>1;
                   %I(I == Imode) = backgroundColor;
                   for dim = 1:size(I,3)
                       temp_I = I(:,:,dim);
                       temp_I(ind) = backgroundColor;
                       I(:,:,dim) = temp_I;
                   end
                   images(:,:,:,ii) = I;
                end   
            end
        end
            
        % Set durations and ISI
        durations = ones(1,size(images,4))*0.5;
        ISI = zeros(1,size(images,4));

        % Generate a number specific for this stimulusType and use
        % this to set seed for stimulus sequence generator below
        % (so we don't use the same sequence for each stimulusType)
        rng('shuffle');
        taskIDs = randi(100);
        taskID = taskIDs(1);
end

% Make individual trial sequences
numberOfStimuli = size(images,4);
% Fix the seed for the random generator such that the same sequence
% will be generated based on the run Number
rng(runNumber+taskID,'twister'); 
stim_seq = randperm(numberOfStimuli);

% Add blank
images(:,:,:,end+1) = mode(images(:));
BLANK = size(images,4);


% This is the stimulus structure used by vistadisp
stimulus              = [];
stimulus.cmap         = stimParams.stimulus.cmap;
stimulus.srcRect      = stimParams.stimulus.srcRect;
stimulus.dstRect      = stimParams.stimulus.destRect;
stimulus.display      = stimParams.display;

% Put everything into stimulus struct
stimulus.categories   = categories;
stimulus.images       = images;
stimulus.im_cell      = im_cell;
stimulus.cat          = catindex;

stimulus.duration     = durations;
stimulus.ISI          = ISI;
stimulus.trialindex   = stim_seq;

% Update durations for temporal stimuli
for ii = 1:numberOfStimuli
    idx = stimulus.trialindex(ii);

    if stimulus.ISI(idx)>0
        stimulus.trial(ii).seqtiming = [...
            [0 stimulus.duration(idx)] ... pulse one
            [0 stimulus.duration(idx)] + stimulus.ISI(idx) + stimulus.duration(idx)... ... pulse two
            ];
        stimulus.trial(ii).seq = [idx BLANK idx BLANK];
    else
        stimulus.trial(ii).seqtiming = [0 stimulus.duration(idx)];
        stimulus.trial(ii).seq = [idx BLANK];
    end
end

% Experiment timing            
fprintf('[%s]: Calculating stimulus timing for: %s\n', mfilename,  site);

% Generate ITIs
numberOfStimuli = size(stimulus.images,4)-1;

switch(lower(stimParams.modality))
    case 'fmri'
        ITI_min  = 3;
        ITI_max  = 6;
        prescan  = round(12/TR)*TR; % seconds
        postscan = prescan; % seconds

        % Jitter ITIs
        ITIs = linspace(ITI_min,ITI_max,numberOfStimuli-1);                

        % Round off to onsetMultiple
        ITIs = round(ITIs/onsetTimeMultiple)*onsetTimeMultiple;

    case {'ecog' 'eeg' 'meg'}
        ITI_min  = 1.25;
        ITI_max  = 1.75;
        prescan  = 3; % seconds
        postscan = 3; % seconds

        % Jitter ITIs
        ITIs = linspace(ITI_min,ITI_max,numberOfStimuli-1);

    otherwise
        error('Unknown modality')
end

stimulus.ITI          = ITIs;
stimulus.prescan      = prescan; % seconds
stimulus.postscan     = postscan; % seconds

% Generate random ITI order
rng('shuffle'); 
iti_seq = randperm(numberOfStimuli-1);

% Compute onsets based on modality-specific ITIs
onsets = cumsum([stimulus.prescan stimulus.ITI(iti_seq)]);

% Match the stimulus presentation to the frame rate
frameRate = stimParams.display.frameRate;
onsets = round(onsets*frameRate)/frameRate;
stimulus.onsets = onsets;

% Put trials together for whole sequence in 'sparse' format: add
% blank at beginning and end, add offsets
seq_sparse       = BLANK; % initialize with blank at time 0
seqtiming_sparse = 0;     % initialize with blank at time 0
for ii = 1:numberOfStimuli
    this_trial_seq = stimulus.trial(ii).seq;
    this_trial_seqtiming = stimulus.trial(ii).seqtiming + onsets(ii);
    seq_sparse = [seq_sparse this_trial_seq];
    seqtiming_sparse = [seqtiming_sparse this_trial_seqtiming];
end    
seq_sparse(end+1)       = BLANK;
seqtiming_sparse(end+1) = seqtiming_sparse(end);

% Put sparse stimulus timing sequences in struct
stimulus.seq_sparse = seq_sparse;
stimulus.seqtiming_sparse = seqtiming_sparse;

% Generate whole sequence at frame Rate resolution 
% Add post-scan stimulus period
%seqtiming = 0:1/frameRate:seqtiming_sparse(end)+max(stimulus.duration)+stimulus.postscan;
seqtiming = 0:1/frameRate:seqtiming_sparse(end)+stimulus.postscan;
seq = zeros(size(seqtiming))+BLANK;
for ii = length(stimulus.seqtiming_sparse):-1:2
    idx = round(seqtiming,4) < round(stimulus.seqtiming_sparse(ii),4);
    seq(idx) = stimulus.seq_sparse(ii-1);
end
seq(end) = stimulus.seq_sparse(end);

% Put interpolated timing sequences in struct
stimulus.seq = seq;
stimulus.seqtiming = seqtiming;

% Add fixation sequence
minDurationInSeconds = 1;
maxDurationInSeconds = 5;
fixSeq = createFixationSequence(stimulus, 1/frameRate, minDurationInSeconds, maxDurationInSeconds);
stimulus.fixSeq = fixSeq;

% Add triggers for non-fMRI modalities
switch lower(stimParams.modality)
    case 'fmri'
        % no triggers for fMRI
    otherwise
        % create an empty trigger sequence
        trigSeq  = zeros(size(stimulus.seq));
        % find the onsets of the stimuli in the sequence
        [~,onsetIndices] = intersect(round(stimulus.seqtiming,4),round(stimulus.onsets,4));
        assert(length(onsetIndices) == length(stimulus.onsets));
        % use the CATEGORICAL labels as trigger codes
        trigSeq(onsetIndices) = stimulus.cat(stimulus.seq(onsetIndices));
        % add task ONSET and OFFSET trigger
        trigSeq(1)   = 256;
        trigSeq(end) = 256;
        stimulus.trigSeq = trigSeq;
end

% Sparsify the stimulus sequence
maxUpdateInterval = 0.25;
stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval);

stimulus.modality = stimParams.modality;

% Generate a save name
fname = sprintf('%s_%s_%d.mat', site, lower(stimulusType), runNumber);

% Add table with elements to write to tsv file for BIDS
onset           = round(stimulus.onsets,3)';
duration        = round(stimulus.duration(stimulus.trialindex),3)';
ISI             = round(stimulus.ISI(stimulus.trialindex),3)';
trial_type      = stimulus.cat(stimulus.trialindex)'; 
trial_name      = stimulus.categories(trial_type - min(stimulus.cat)+1)';
stim_file       = repmat(fname, numberOfStimuli ,1);
stim_file_index = stimulus.trialindex';

stimulus.tsv    = table(onset, duration, ISI, trial_type, trial_name, stim_file, stim_file_index);

stimulus.site   = site;

% save 
fprintf('[%s]: Saving stimuli in: %s\n', mfilename, fullfile(vistadispRootPath, 'StimFiles',  fname));
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus', '-v7.3')

return


