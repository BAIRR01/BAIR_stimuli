function stimMakeSpatiotemporalExperiment(stimParams, runID, stimulusType)

%% SPATIO TEMPORAL
% For visual experiments, we use band-pass, gray-scale images, spanning
% many stimulus dimensions. Twelve were used in a prior publication [69,70],
% varying in contrast, number of component orientations (1, 2 or 16
% superimposed gratings), or spacing between contrast elements (from very
% sparse to very dense). Twelve are natural images of faces, objects, and
% scenes (also gray-scale, band-pass). These stimuli will be presented for
% 0.5 seconds each. Twelve other stimuli are simple noise patterns shown
% with different temporal profiles (single pulses with variable duration;
% or multiple pulses with variable interstimulus interval).

% Runs are divided up into 3 STIMULUS TYPES:
    % Spatial : Patterns (noise patterns varying in contrast and density; gratings) 
    % Spatial : Objects  (letters, faces, houses)
    % Temporal: Patterns (differing in duration and ISI)
    
% 24 unique exemplars divided across two unique RUNS for each TYPE
% Per RUN we have 36 stimulus presentations in total:

% CRF      - noise patterns at 5 contrast levels (5*3)                     
% Orient   - gratings at 3 complexity levels (single grating, plaid, circular) (3*3)
% Sparsity - noise patterns at 4 density levels (5th level included in CRF) (4*3)

% Faces    - faces regenerated from database used for Kay et al., 2013 (1*12)
% Letters  - letters regenerated from database used for Kay et al., 2013 (1*12)
% Houses   - houses regenerated from database used for Kay et al., 2013 (1*12)

% 1 Pulse  - noise patterns at 6 durations each (6*3)
% 2 Pulses - noise patterns at 6 ISIs (fixed duration) (6*3)

% Order of exemplars fixed between EcoG, fMRI, MEG.

%% Make the images

% Determine if we're creating the master or loading & resizing for a specific display
site = stimParams.experimentSpecs.sites{1};
imageSizeInPixels = size(stimParams.stimulus.images);

% This is the stimulus structure used by vistadisp
stimulus.cmap         = stimParams.stimulus.cmap;
stimulus.srcRect      = stimParams.stimulus.srcRect;
stimulus.dstRect      = stimParams.stimulus.destRect;
stimulus.display      = stimParams.display;

switch site
    case 'Master'
        
        switch stimulusType
            case 'SPATIALPATTERN'
                
                fprintf('[%s]: Creating Master stimuli for stimulusType: SPATIALPATTERN, runID: %d.\n',mfilename, runID);

                categories = {...
                    'CRF-1' ...
                    'CRF-2' ...
                    'CRF-3' ...
                    'CRF-4' ...
                    'CRF-5' ...
                    'GRATING' ...
                    'PLAID' ...
                    'CIRCULAR' ...
                    'SPARSITY-1' ...
                    'SPARSITY-2' ...
                    'SPARSITY-3' ...
                    'SPARSITY-4' ...
                    };
            
                numberOfImagesPerCat = 3;
                
                % Pre-allocate arrays to store images
                images = uint8(zeros([imageSizeInPixels length(categories) * numberOfImagesPerCat]));
                im_cell = cell([1 length(categories)]);
                imCount = 1;

                % Create CRF stimuli
                fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: CRF.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

                % From Kay et al, 2013 PLOS CB
                %   These stimuli were constructed by varying the contrast
                %   of the noise patterns used in SPACE. Ten different
                %   contrast levels were used: 1%, 2%, 3%, 4%, 6%, 9%, 14%,
                %   21%, 32%, and 50%. These contrast levels are relative
                %   to the contrast of the patterns used in SPACE, which is
                %   taken to be 100%.

                % Category-specific settings
                categoryIndex = find(contains(categories, {'CRF'}));
                numberOfCategories = length(categoryIndex);
                
                contrastLevels = [0.0625 0.125 0.25 0.5 1]; % or should we use log/linear?
                densityLevel = 20; % middle density stimulus for SPARSITY

                % Create the stimuli
                for cc = 1:numberOfCategories
                    for ii = 1:numberOfImagesPerCat
                        imageForThisTrial = createPatternStimulus(stimParams, densityLevel, contrastLevels(cc));

                        % Double to unsigned 8 bit integer, needed for vistadisp
                        % Note: luminance values need to range between -0.5 and 0.5
                        % prior to 8Bit conversion to avoid clipping
                        image8Bit = uint8((imageForThisTrial+.5)*255);
                        images(:,:,imCount) = image8Bit;
                        im_cell{categoryIndex(cc)}(:,:,ii) = image8Bit;
                        imCount = imCount + 1;
                    end
                end

                % Create GRATING / PLAID / CIRCULAR
                fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: GRATING PLAID CIRCULAR.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

                % From Kay et al, 2013 PLOS CB
                % GRATING (4 stimuli). These stimuli consisted of horizontal
                %   sinusoidal gratings at 2%, 4%, 9%, and 20% Michelson contrast.
                %   The spatial frequency of the gratings was fixed at 3 cycles per
                %   degree. Each stimulus consisted of gratings with the same
                %   contrast but nine different phases (equally spaced from 0 to 2p).
                % PLAID (4 stimuli). These stimuli consisted of plaids at 2%, 4%,
                %   9%, and 20% contrast (defined below). Each condition comprised
                %   nine plaids, and each plaid was constructed as the sum of a
                %   horizontal and a vertical sinusoidal grating (spatial frequency 3
                %   cycles per degree, random phase). The plaids were scaled in
                %   contrast to match the root-mean-square (RMS) contrast of the
                %   GRATING stimuli. For example, the plaids in the 9% condition were
                %   scaled such that the average RMS contrast of the plaids is
                %   identical to the average RMS contrast of the gratings in the 9%
                %   GRATING stimulus. 
                % CIRCULAR (4 stimuli). These stimuli were identical to the
                %   PLAID stimuli except that sixteen different orientations
                %   were used instead of two.

                % Category-specific settings
                categoryIndex = find(contains(categories,  {'GRATING', 'PLAID', 'CIRCULAR'}));
                numberOfCategories = length(categoryIndex);
            
                cyclesPerDegree     = 3;     % cycles per degree
                peak2peakContrast   = 50;    % percentage Michelson contrast

                % Create the stimuli
                for cc = 1:numberOfCategories
                    switch categories{categoryIndex(cc)}
                        case 'GRATING'
                            gratingOrientation = pi/2; % horizontal
                        case 'PLAID'
                            gratingOrientation = [pi/2 pi]; % horizontal + vertical
                        case 'CIRCULAR'
                            gratingOrientation =  (1:16)/pi; % 16 superimposed
                    end
                    for ii = 1:numberOfImagesPerCat
                        imageForThisTrial = createGratingStimulus(stimParams,cyclesPerDegree,gratingOrientation, peak2peakContrast);

                        % Double to unsigned 8 bit integer, needed for vistadisp
                        image8Bit = uint8((imageForThisTrial+.5)*255);
                        images(:,:,imCount) = image8Bit;
                        im_cell{categoryIndex(cc)}(:,:,ii) = image8Bit;
                        imCount = imCount + 1;
                    end
                end

                 % Create SPARSITY 
                fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: SPARSITY.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

                % From Kay et al, 2013 PLOS CB
                %   We generated noise patterns using cutoff frequencies of
                %   2.8, 1.6, 0.9, 0.5, and 0.3 cycles per degree, and
                %   numbered these from 1 (smallest separation) to 5
                %   (largest separation). The noise patterns used in SPACE
                %   correspond to separation 4; thus, we only constructed
                %   stimuli for the remaining separations 1, 2, 3, and 5.
                %   The noise patterns occupied the full stimulus extent
                %   (no aperture masking).

                % Category-specific settings
                categoryIndex = find(contains(categories, {'SPARSITY'}));
                numberOfCategories = length(categoryIndex);
            
                contrastLevel = 1; % fixed to max
                densityLevels = [10 15 25 30]; % middle density is CRF full contrast stim

                % Create the stimuli
                for cc = 1:numberOfCategories
                    for ii = 1:numberOfImagesPerCat
                        imageForThisTrial = createPatternStimulus(stimParams, densityLevels(cc), contrastLevel);

                        % Double to unsigned 8 bit integer, needed for vistadisp
                        image8Bit = uint8((imageForThisTrial+.5)*255);
                        images(:,:,imCount) = image8Bit;
                        im_cell{categoryIndex(cc)}(:,:,ii) = image8Bit;
                        imCount = imCount + 1;
                    end
                end
                
                % Set durations and ISI
                durations = ones(1,size(images,3))*0.5;
                ISI = zeros(1,size(images,3));
                
                % Generate a number specific for this stimulusType and use
                % this to set seed for stimulus sequency generator below
                % (so we don't use the same sequence for each stimulusType)
                taskID = 10; 
                
            case 'SPATIALOBJECT'
                
                categories = {...
                    'FACES' ...
                    'LETTERS' ...
                    'HOUSES' ...
                    };
            
                numberOfImagesPerCat = 12;

                % Pre-allocate arrays to store images
                images = uint8(zeros([imageSizeInPixels length(categories) * numberOfImagesPerCat]));
                im_cell = cell([1 length(categories)]);
                imCount = 1;

                % Create FACES / LETTERS / HOUSES
                fprintf('[%s]: Creating Master stimuli for stimulusType: SPATIALOBJECT.\n',mfilename);

                % From Kay et al, 2013 PLOS CB
                %   Each image was whitened (to remove low-frequency bias)
                %   and then filtered with the custom band-pass filter
                %   (described previously). Finally, the contrast of each
                %   image was scaled to fill the full luminance range.

                % Load original, unfiltered face, house, letter stimuli 
                % CHANGE TO DOWNLOAD FROM WIKI??
                
                % Download spatiotemporal stimuli
                %         images = []; for ii = 1:13 % why 13?
                %             readPth =
                %             sprintf('https://wikis.nyu.edu/download/attachments/85394548/spatiotemporal%d.mat?api=v2',
                %             ii); stimDir = fullfile(BAIRRootPath,
                %             'stimuli'); fname =
                %             sprintf('spatiotemporal%d.mat', ii); writePth
                %             = fullfile(stimDir, fname); if
                %             ~exist(writePth, 'file'),
                %             websave(writePth,readPth); end tmp =
                %             load(writePth); images = [images tmp.im];
                %         end

                fprintf('[%s]: Loading original, unfiltered stimuli...\n',mfilename);
                load('/Volumes/server/Projects/BAIR/Stimuli/Kay2013_unfilteredstimuli/stim.mat');

                % Category-specific settings
                categoryIndex = find(contains(categories,  {'FACES', 'LETTERS', 'HOUSES'}));
                numberOfCategories = length(categoryIndex);
                
                imageProcessingParams.imageScaleFactor  = 0.4; % size of FHL stimuli relative to the pattern stimuli
                imageProcessingParams.maskScaleFactor   = 0.85; % size of mask relative to filtered image
                imageProcessingParams.maskSoftEdge      = 1.1; % size of 'second radius' of mask (determines smoothness of mask edge)
                imageProcessingParams.contrastRange     = 1.5; % range to scale the contrast prior to 8Bit conversion (a range of 1 results in no clipping)
                
                % Create the stimuli
                for cc = 1:numberOfCategories
                    
                    switch categories{categoryIndex(cc)}
                        case 'FACES'
                            imageArray = faces;
                            fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: FACES.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));
                        case 'LETTERS'
                            imageArray = letters;
                            fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: LETTERS.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));
                            % For LETTERS ONLY: Crop and resize letter
                            % array prior to filtering to get thicker lines
                            origSizeinPixels = size(letters);
                            cropOrigin = [(origSizeinPixels(1)+1)./2 (origSizeinPixels(1)+1)./2];
                            cropRadius = [-0.3*origSizeinPixels(1)/2 0.3*origSizeinPixels(1)/2]; % optimal crop radius determined empirically REDO
                            cropIndex  = round(cropOrigin+cropRadius);
                            imageArray = imageArray(cropIndex(1):cropIndex(2), cropIndex(1):cropIndex(2),:);
                            imageArray = imresize(imageArray,origSizeinPixels([1 2]));               
                        case 'HOUSES'
                            imageArray = houses;                       
                            fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: HOUSES.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));
                    end
                    
                    % Pick which stimuli to select from original set
                    % ODD for runNum == 1, EVEN for runNum == 2;
                    imageIndex = runID:2:numberOfImagesPerCat*2+runID-1;
                    % Should we balance gender/nationalities for the faces?
                    % Should we remove certain faces (eg Jon?)
               
                    for ii = 1:numberOfImagesPerCat
                        inputImage = imageArray(:,:,imageIndex(ii));
                        imageForThisTrial = createFilteredStimulus(stimParams,inputImage,imageProcessingParams);

                        % Double to unsigned 8 bit integer, needed for vistadisp
                        image8Bit = uint8((imageForThisTrial+.5)*255);
                        images(:,:,imCount) = image8Bit;
                        im_cell{categoryIndex(cc)}(:,:,ii) = image8Bit;
                        imCount = imCount + 1;
                    end
                end

                % Set durations and ISI
                durations = ones(1,size(images,3))*0.5;
                ISI = zeros(1,size(images,3));
            
                % Generate a number specific for this stimulusType and use
                % this to set seed for stimulus sequency generator below
                % (so we don't use the same sequence for each stimulusType)
                taskID = 11; 
                
            case 'TEMPORALPATTERN'
                 
                fprintf('[%s]: Creating Master stimuli for stimulusType: TEMPORALPATTERN.\n',mfilename);

                categories = {...
                    'ONEPULSE-1' ...
                    'ONEPULSE-2' ...
                    'ONEPULSE-3' ...
                    'ONEPULSE-4' ...
                    'ONEPULSE-5' ...
                    'ONEPULSE-6' ...
                    'TWOPULSE-1' ...
                    'TWOPULSE-2' ...
                    'TWOPULSE-3' ...
                    'TWOPULSE-4' ...
                    'TWOPULSE-5' ...
                    'TWOPULSE-6' ...
                    };
                
                numberOfImagesPerCat = 3;
                
                % Pre-allocate arrays to store images
                images = uint8(zeros([imageSizeInPixels length(categories) * numberOfImagesPerCat]));
                im_cell = cell([1 length(categories)]);
                imCount = 1;

                % Create TEMPORAL stimuli
                fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: ONE PULSE.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

                % Category-specific settings
                numberOfCategories = length(categories);
                categoryIndex = 1:length(categories);
                
                contrastLevel = 1; % max contrast
                densityLevel = 20; % middle density

                % Create the stimuli
                for cc = 1:numberOfCategories
                    if cc == (numberOfCategories/2)+1
                        fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: TWO PULSE.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));
                    end
                    for ii = 1:numberOfImagesPerCat
                        imageForThisTrial = createPatternStimulus(stimParams, densityLevel, contrastLevel);

                        % Double to unsigned 8 bit integer, needed for vistadisp
                        image8Bit = uint8((imageForThisTrial+.5)*255);
                        images(:,:,imCount) = image8Bit;
                        im_cell{categoryIndex(cc)}(:,:,ii) = image8Bit;
                        imCount = imCount + 1;
                    end
                end
                
                % Set durations and ISIs           
                tempIndex = [1 2 4 8 16 32]/60;
                
                durations = [];
                % One pulse durations:
                for ii = 1:length(tempIndex)
                    durations = [durations ones(1,numberOfImagesPerCat)*tempIndex(ii)];
                end
                
                % Append two pulse durations:
                durations = [durations ones(1,length(tempIndex)*numberOfImagesPerCat)*8/60];
                
                % One pulse ISI:
                ISI = zeros(1,18);

                % Append two pulse ISI:
                for ii = 1:length(tempIndex)
                    ISI = [ISI ones(1,numberOfImagesPerCat)*tempIndex(ii)];
                end    
                
                % Generate a number specific for this stimulusType and use
                % this to set seed for stimulus sequency generator below
                % (so we don't use the same sequence for each stimulusType)
                taskID = 12; 
        end
        
        % DEBUG: plot generated images
        figure;hold on
        for ii = 1:size(images,3)
            subplot(ceil(sqrt(size(images,3))),ceil(sqrt(size(images,3))),ii);
            imshow(images(:,:,ii));
        end

        % Make individual trial sequences
        numberOfStimuli = size(images,3);
        % Fix the seed for the random generator such that the same sequence
        % will be generated based on the run Number
        rng(runID+taskID,'twister'); 
        stim_seq = randperm(numberOfStimuli);
        
        % Add blank
        images(:,:,end+1) = mode(images(:));
        BLANK = size(images,3);

        % Put everything into stimulus struct
        stimulus.categories   = categories;
        stimulus.im_cell      = im_cell;
        stimulus.images       = images;

        stimulus.duration     = durations;
        stimulus.ISI          = ISI;
        stimulus.cat          = stim_seq;
        
        % Update durations for temporal stimuli
        for ii = 1:numberOfStimuli
            idx = stim_seq(ii);
        
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
      
        otherwise    

            % Resize the Master stimuli to the required stimulus size for this
                    % modality and display
            fprintf('[%s]: Loading Master stimuli for stimulusType: %s, runID: %d \n',mfilename, stimulusType, runID);

            % Load the Master stimuli
            stimulus = loadBAIRStimulus(stimulusType, 'Master', runID);

            % Resize         
            fprintf('[%s]: Resizing Master stimuli for: %s\n', mfilename,  site);
            imageSizeInPixels = size(stimParams.stimulus.images);
            images = imresize(stimulus.images, imageSizeInPixels);

            im_cell = cell(size(stimulus.im_cell));
            for ii = 1:size(stimulus.im_cell,2)
                im_cell{ii} = imresize(stimulus.im_cell{ii}, imageSizeInPixels);
            end

            % Soft circular mask (1 pixel of blurring per 250 pixels in the image)
            supportDiameter       = imageSizeInPixels(1);
            maskRadius            = (stimParams.stimulus.srcRect(3) - stimParams.stimulus.srcRect(1))/2;
            circularMask          = mkDisc(supportDiameter, maskRadius, (imageSizeInPixels+1)./2, 1/250 * imageSizeInPixels(1));
            imagesDouble          = double(images)/255-.5;
            imagesMasked          = bsxfun(@times,imagesDouble,circularMask);
            images                = uint8((imagesMasked+.5)*255);

            % Overwrite Master stimuli with resized images
            stimulus.im_cell      = im_cell;
            stimulus.images       = images;   
            
            % DEBUG: plot generated images
            figure;hold on
            for ii = 1:size(images,3)-1
                subplot(ceil(sqrt(size(images,3)-1)),ceil(sqrt(size(images,3)-1)),ii);
                imshow(images(:,:,ii));
            end
            
            % Experiment timing 
            
            fprintf('[%s]: Calculating stimulus timing for: %s\n', mfilename,  site);

            switch(lower(stimParams.modality))
                case 'fmri'
                    ITI_min  = 3;
                    ITI_max  = 6;
                    prescan  = 9; % seconds
                    postscan = 9; % seconds
                case {'ecog' 'eeg' 'meg'}
                    ITI_min  = 1.25;
                    ITI_max  = 1.75;
                    prescan  = 2; % seconds
                    postscan = 2; % seconds
                otherwise
                    error('Unknown modality')
            end

            % Generate ITIs
            numberOfStimuli = size(stimulus.images,3)-1;
            ITIs = linspace(ITI_min,ITI_max,numberOfStimuli);

            stimulus.ITI          = ITIs;
            stimulus.prescan      = prescan; % seconds
            stimulus.postscan     = postscan; % seconds

            % Generate random ITI order
            rng('shuffle'); 
            iti_seq = randperm(numberOfStimuli);

            % Compute onsets based on modality-specific ITIs
            stimulus.onsets = cumsum([stimulus.prescan stimulus.ITI(iti_seq)]);
            stimulus.onsets = stimulus.onsets(1:end-1);

            BLANK = size(stimulus.images,3);
            stimulus.seq       = BLANK; % initialize with blank at time 0
            stimulus.seqtiming = 0;     % initialize with blank at time 0

            % Put trials together for whole sequence
            trigSeq = 0; % initialize trigger sequence with 0
            for ii = 1:numberOfStimuli
                this_trial_seq = stimulus.trial(ii).seq;
                this_trial_seqtiming = stimulus.trial(ii).seqtiming + stimulus.onsets(ii);
                stimulus.seq = [stimulus.seq this_trial_seq];
                stimulus.seqtiming = [stimulus.seqtiming this_trial_seqtiming];

                this_trial_trig_seq = zeros(size(this_trial_seq));
                this_trial_trig_seq(1) = 1;
                trigSeq   = [trigSeq this_trial_trig_seq];
            end
            stimulus.seq(end+1) = BLANK;
            stimulus.seqtiming(end+1) = stimulus.seqtiming(end) + stimulus.postscan;

            % Interpolate to 60 frames / second
            seqtiming = 0:1/60:stimulus.seqtiming(end);
            seq = zeros(size(seqtiming));

            for ii = length(stimulus.seqtiming):-1:2
                idx = seqtiming < stimulus.seqtiming(ii);
                seq(idx) = stimulus.seq(ii-1);
            end
            seq(end) = stimulus.seq(end);

            % Put full and sparse timing sequences in struct
            stimulus.seqtiming_sparse = stimulus.seqtiming;
            stimulus.seq_sparse = stimulus.seq;
            stimulus.seq = seq;
            stimulus.seqtiming = seqtiming;

            % Triggers
            trigSeq  = zeros(size(stimulus.seq));
            diodeSeq = zeros(size(stimulus.seq)); 

            for ii = 1:length(stimulus.onsets)
                [~, idx] = min(abs(stimulus.seqtiming-stimulus.onsets(ii)));
                trigSeq(idx) = stimulus.cat(ii);
                diodeSeq(idx) = 1;
            end

            switch lower(stimParams.modality)
                case 'fmri'
                otherwise
                    stimulus.trigSeq = trigSeq;
                    stimulus.diodeSeq = diodeSeq;
            end
   
            stimulus.modality = stimParams.modality;
            stimulus.site     = site;
end

% save 
fname = sprintf('%s_%s_%d.mat', lower(stimulusType), site, runID);
fprintf('[%s]: Saving Master stimuli as: %s\n', mfilename, fname);

save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus', '-v7.3')

return


% OLD:
tmp = cumsum(cellfun(@length, whichIm));
whichIm_Idx = [[0 tmp(1:end-1)]+1; tmp];

stim_seq = randperm(num_cats);
for ii = 1:num_cats
    idx = stim_seq(ii);

    % choose exemplar based on run number
    possible_exemplars = whichIm_Idx(1,idx):whichIm_Idx(2,idx);
    n = length(possible_exemplars);
    this_exemplar = mod(runNum, n);

    thisim = whichIm_Idx(1,idx)+this_exemplar;

    if stimulus.ISI(idx)>0
        stimulus.trial(ii).seqtiming = [...
            [0 stimulus.duration(idx)] ... pulse one
            [0 stimulus.duration(idx)] + stimulus.ISI(idx) + stimulus.duration(idx)... ... pulse two
            ];
        stimulus.trial(ii).seq = [thisim BLANK thisim BLANK];
    else
        stimulus.trial(ii).seqtiming = [0 stimulus.duration(idx)];
        stimulus.trial(ii).seq = [thisim BLANK];
    end

end