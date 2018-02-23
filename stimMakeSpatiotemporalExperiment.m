function stimMakeSpatiotemporalExperiment(stimParams, runNum)

%% SPATIO TEMPORAL (12 repeats, ECoG, fMRI; 24 for E/MEG !?)
% For visual experiments, we use band-pass, gray-scale images, spanning
% many stimulus dimensions. Twelve were used in a prior publication [69,70],
% varying in contrast, number of component orientations (1, 2 or 16
% superimposed gratings), or spacing between contrast elements (from very
% sparse to very dense). Twelve are natural images of faces, objects, and
% scenes (also gray-scale, band-pass). These stimuli will be presented for
% 0.5 seconds each. Twelve other stimuli are simple noise patterns shown
% with different temporal profiles (single pulses with variable duration;
% or multiple pulses with variable interstimulus interval).

% CRF      - 5 (zebra)                     KNK 162 164 166 168 116
% Orient   - 3 (grating, plaid, circular)  KNK 150, 154, 158 (*HC)
% Sparsity - 4 (zebras)                    KNK 181 182 183 184
% 1 Pulse  - 6 (zebra??)                   KNK 183 * 6
% 2 Pulses - 6 (zebra??)                   KNK 183 * 6

% Faces -    4                             KNK 171 (sample 6 * 8 for 12 runs, 4 each)
% Letters -  4                             KNK 173 (sample 6 * 8 for 12 runs, 4 each)
% Scenes -   4                             KNK 175 (sample 6 * 8 for 12 runs, 4 each)

%% Make the images

% Determine if we're creating the master or loading & resizing for a specific display
site = stimParams.experimentSpecs.sites{1};
imageSizeInPixels = size(stimParams.stimulus.images);

% This is the stimulus structure used by vistadisp
stimulus.cmap         = stimParams.stimulus.cmap;
stimulus.srcRect      = stimParams.stimulus.srcRect;
stimulus.dstRect      = stimParams.stimulus.destRect;

switch site
    case 'Master'
        
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
            'FACES-1' ...
            'FACES-2' ...
            'FACES-3' ...
            'FACES-4' ...
            'LETTERS-1' ...
            'LETTERS-2' ...
            'LETTERS-3' ...
            'LETTERS-4' ...
            'HOUSES-1' ...
            'HOUSES-2' ...
            'HOUSES-3' ...
            'HOUSES-4' ...
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
        
        % Create CRF stimuli
        fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: CRF.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

        % From Kay et al, 2013 PLOS CB
%         These stimuli were constructed by varying the contrast of the
%         noise patterns used in SPACE. Ten different contrast levels were
%         used: 1%, 2%, 3%, 4%, 6%, 9%, 14%, 21%, 32%, and 50%. These
%         contrast levels are relative to the contrast of the patterns used
%         in SPACE, which is taken to be 100%.
        
        % Preallocate arrays to store stimuli (note: we can do this more
        % efficiently by preallocating for the entire set at the beginning
        % and fill it with the various categories (e.g. CRF, grating), as
        % we go along. For now, I like to be able to keep track of those
        % separately in category-specific arrays, and then concatenate the
        % category-specific arrays at the end.
        numberOfCategories = length(find(contains(categories, 'CRF')));
        numberOfImagesPerCat = 6;
        CRFimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        CRFim_cell = cell([1 numberOfCategories]);
        
        % Category-specific settings
        contrastLevels = [0.0625 0.125 0.25 0.5 1]; % or should we use log/linear?
        densityLevels = 20; % middle density stimulus for SPARSITY

        % Create the stimuli
        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createPatternStimulus(stimParams, densityLevels, contrastLevels(cc));
                
                % Double to unsigned 8 bit integer, needed for vistadisp
                % Note: luminance values need to range between -0.5 and 0.5
                % prior to 8Bit conversion to avoid clipping
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                CRFimages(:,:,imCount) = image8Bit;
                
            end
            CRFim_cell{cc} = CRFimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

        % Create SPARSITY 
        fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: SPARSITY.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

        % From Kay et al, 2013 PLOS CB
%         We generated noise patterns using cutoff frequencies of 2.8, 1.6,
%         0.9, 0.5, and 0.3 cycles per degree, and numbered these from 1
%         (smallest separation) to 5 (largest separation). The noise
%         patterns used in SPACE correspond to separation 4; thus, we only
%         constructed stimuli for the remaining separations 1, 2, 3, and 5.
%         The noise patterns occupied the full stimulus extent (no aperture
%         masking).

        % Preallocate arrays to store stimuli
        numberOfCategories = length(find(contains(categories, 'SPARSITY')));
        numberOfImagesPerCat = 6;
        SPARSITYimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        SPARSITYim_cell = cell([1 numberOfCategories]);

        % Category-specific settings
        contrastLevels = 1; % or should we use log/linear?
        densityLevels = [10 15 25 30]; % middle density is CRF full contrast stim

        % Create the stimuli
        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createPatternStimulus(stimParams, densityLevels(cc), contrastLevels);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                SPARSITYimages(:,:,imCount) = image8Bit;

            end
            SPARSITYim_cell{cc} = SPARSITYimages(:,:,(cc-1)*numberOfImagesPerCat+1);
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
        
        % Preallocate arrays to store stimuli
        numberOfCategories = length(find(contains(categories, {'GRATING', 'PLAID', 'CIRCULAR'})));
        numberOfImagesPerCat = 6;
        GPCimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        GPCim_cell = cell([1 numberOfCategories]);

        % Category-specific settings
        cyclesPerDegree     = 3;     % cycles per degree
        peak2peakContrast   = 50;    % percentage Michelson contrast
        
        % Create the stimuli
        for cc = 1:numberOfCategories
            if cc == 1
                gratingOrientation = pi/2; % GRATING: horizontal
            elseif cc == 2
                gratingOrientation = [pi/2 pi]; % PLAID: horizontal + vertical
            elseif cc == 3
                gratingOrientation =  (1:16)/pi; % CIRCULAR: 16 superimposed
            end
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createGratingStimulus(stimParams,cyclesPerDegree,gratingOrientation, peak2peakContrast);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                GPCimages(:,:,imCount) = image8Bit;

            end
            GPCim_cell{cc} = GPCimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

        % DEBUG  % compare with old version
%         load('/Users/winawerlab/Box Sync/Stimuli/spatiotemporal_NYU-ECOG_101.mat')

%         figure;imshow(stimulus.images(:,:,31)); title('OLD');
%         figure;imshow(stimulus.images(:,:,37)); title('OLD');
%         figure;imshow(stimulus.images(:,:,43)); title('OLD');

        % Create FACES / LETTERS / HOUSES
        fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: FACES LETTERS HOUSES.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

        % From Kay et al, 2013 PLOS CB
%          Each image was whitened (to remove low-frequency bias) and then
%          filtered with the custom band-pass filter (described
%          previously). Finally, the contrast of each image was scaled to
%          fill the full luminance range.

        % Preallocate arrays to store stimuli
        numberOfCategories = length(find(contains(categories, {'FACES', 'LETTERS', 'HOUSES'})));
        numberOfImagesPerCat = 2;
        FLHimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        FLHim_cell = cell([1 numberOfCategories]);
     
        % Load original, unfiltered face, house, letter stimuli 
        % CHANGE TO DOWNLOAD FROM WIKI
        load('/Volumes/server/Projects/BAIR/Stimuli/Kay2013_unfilteredstimuli/stim.mat');
        
        % Category-specific settings
        imageProcessingParams.imageScaleFactor  = 0.4; % size of FHL stimuli relative to the pattern stimuli
        imageProcessingParams.maskScaleFactor   = 0.85; % size of mask relative to filtered image
        imageProcessingParams.maskSoftEdge      = 1.1; % size of 'second radius' of mask (determines smoothness of mask edge)
        imageProcessingParams.contrastRange     = 1.5; % range to scale the contrast prior to 8Bit conversion (a range of 1 results in no clipping)
        imageIndex = [30 35; % FACES
                      24 83; 
                      74 13; 
                      59 41;
                      13 21; % LETTERS
                      6 20; 
                      1 9; 
                      8 24;
                      33 7;  % HOUSES
                      24 20; 
                      41 3; 
                      31 25]; 

        % Create the stimuli
        for cc = 1:numberOfCategories    
            
            % Select images from the right array
            if ismember(cc,1:4) % FACES              
                imageArray = faces;
            elseif ismember(cc,5:8) % LETTERS
                imageArray = letters;
                % For LETTERS: Crop and resize letter array prior to
                % filtering to get thicker lines
                origSizeinPixels = size(letters);
                cropOrigin = [(origSizeinPixels(1)+1)./2 (origSizeinPixels(1)+1)./2];
                cropRadius = [-0.25*origSizeinPixels(1)/2 0.25*origSizeinPixels(1)/2]; % optimal crop radius determined empirically
                cropIndex  = round(cropOrigin+cropRadius);
                imageArray = imageArray(cropIndex(1):cropIndex(2), cropIndex(1):cropIndex(2),:);
                imageArray = imresize(imageArray,origSizeinPixels([1 2]));               
            elseif ismember(cc,9:12) % HOUSES
                imageArray = houses;           
            end
            
            for ii = 1:numberOfImagesPerCat
                inputImage = imageArray(:,:,imageIndex(cc,ii));
                imageForThisTrial = createFilteredStimulus(stimParams,inputImage,imageProcessingParams);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                FLHimages(:,:,imCount) = image8Bit;

            end
            FLHim_cell{cc} = FLHimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end
        
%         % DEBUG  % compare with old version
%         load('/Users/winawerlab/Box Sync/Stimuli/spatiotemporal_NYU-ECOG_101.mat')
%         figure;hold on;
%         pli = 1;
%         for ii = 73:96
%             subplot(3,8,pli);
%             imshow(stimulus.images(:,:,ii),[])
%             pli = pli+1;
%         end

        % Create TEMPORAL
        fprintf('[%s]: Creating Master stimuli at %d x %d pixels resolution: TEMPORAL.\n',mfilename,imageSizeInPixels(1),imageSizeInPixels(2));

        % First and second for double pulse should be same stimulus (TO DO
        % check how this is determined in stim sequencing)
        
        % Preallocate arrays to store stimuli
        numberOfCategories = length(find(contains(categories, 'PULSE')));
        numberOfImagesPerCat = 6;
        PULSEimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        PULSEim_cell = cell([1 numberOfCategories]);
        
        % Category-specific settings
        contrastLevels = 1; % max contrast
        densityLevels = 20; % middle density
        
        % Create the stimuli
        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createPatternStimulus(stimParams, densityLevels, contrastLevels);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                PULSEimages(:,:,imCount) = image8Bit;

            end
            PULSEim_cell{cc} = PULSEimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

        % Concatenate all categories
        images = cat(3, CRFimages, GPCimages, SPARSITYimages, FLHimages, PULSEimages);
        im_cell = cat(2, CRFim_cell, GPCim_cell, SPARSITYim_cell, FLHim_cell, PULSEim_cell);
        
        % Set durations
        durations = [ ...
            ones(1,24)*0.5      ... spatial
            [1 2 4 8 16 32]/60  ... one pulse
            ones(1,6)*8/60      ... two pulse
            ];

        ISI = [ ...
            zeros(1,24)         ... spatial
            zeros(1,6)          ... one pulse
            [1 2 4 8 16 32]/60  ... two pulse
            ];
        
        % Add blank
        images(:,:,end+1) = mode(images(:));
        BLANK = size(images,3);
        
        % Put everything into stimulus struct
        stimulus.categories   = categories;
        stimulus.im_cell      = im_cell;
        stimulus.images       = images;

        stimulus.duration     = durations;
        stimulus.ISI          = ISI;
        
    otherwise    
        
        % Resize the Master stimuli to the required stimulus size for this
                % modality and display
        fprintf('[%s]: Loading Master stimuli for: %s\n', mfilename, site);

        % Load the Master stimuli
        stimulus = loadBAIRStimulus('spatiotemporal', 'Master', runNum);
        
        % Resize         
        fprintf('[%s]: Resizing Master stimuli for: %s\n', mfilename,  site);
        imageSizeInPixels = size(stimParams.stimulus.images);
        images = imresize(stimulus.images, imageSizeInPixels );
        
        im_cell = cell(num_cats);
        for ii = 1:num_cats
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
end

% Experiment timing 
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

num_cats = length(stimulus.categories);
ITIs = linspace(ITI_min,ITI_max,num_cats);

stimulus.ITI          = ITIs;
stimulus.prescan      = prescan; % seconds
stimulus.postscan     = postscan; % seconds

% TO DO: make individual trial sequences
%   * randomize image exemplars across runs? (cf Dora comments)
%   * chunk categories within separate runs for fMRI?

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

% Put trials together for whole sequence
%   randomize ITI order
iti_seq = randperm(num_cats);

stimulus.onsets = cumsum([stimulus.prescan stimulus.ITI(iti_seq)]);
stimulus.onsets = stimulus.onsets(1:end-1);

stimulus.seq       = BLANK; % initialize with blank at time 0
stimulus.seqtiming = 0;     % initialize with blank at time 0

trigSeq   = 0; % initialize trigger sequence with 0

for ii = 1:num_cats
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

seqtiming =  0:1/60:stimulus.seqtiming(end);
seq = zeros(size(seqtiming));

for ii = length(stimulus.seqtiming):-1:2
    idx = seqtiming < stimulus.seqtiming(ii);
    seq(idx) = stimulus.seq(ii-1);
end
seq(end) = stimulus.seq(end);


stimulus.seqtiming_sparse = stimulus.seqtiming;
stimulus.seq_sparse = stimulus.seq;
stimulus.seq = seq;
stimulus.seqtiming = seqtiming;

% triggers
trigSeq  = zeros(size(stimulus.seq));
diodeSeq = zeros(size(stimulus.seq)); % don't need this anymore?

for ii = 1:length(stimulus.onsets)
    [~, idx] = min(abs(stimulus.seqtiming-stimulus.onsets(ii)));
    trigSeq(idx) = stim_seq(ii);
    diodeSeq(idx) = 1;
end

stimulus.cat = stim_seq;
switch lower(stimParams.modality)
    case 'fmri'
    otherwise
        stimulus.trigSeq = trigSeq;
        stimulus.diodeSeq = diodeSeq;
end

stimulus.display  = stimParams.display;
stimulus.modality = stimParams.modality;
stimulus.site     = site;

% save 
fname = sprintf('spatiotemporal_%s_%d', site, runNum);
save(fullfile(vistadispRootPath, 'StimFiles',  fname), 'stimulus', '-v7.3')

return

%% DEBUG 

%% Plot individual stimuli for a specific category ('images')
whichCat = 'FLH'; %pick: CRF, GPC, SPARSITY, FLH, or PULSE;

images = eval([whichCat 'images']);

figure;hold on
for ii = 1:size(images,3)
    subplot(ceil(sqrt(size(images,3))),ceil(sqrt(size(images,3))),ii);
    imshow(images(:,:,ii));
end

%% Plot exemplars for all categories ('im_cell')

figure('Name', 'NEW STIM');hold on
for ii = 1:24%length(stimulus.im_cell)
    subplot(4,6,ii);
    imshow(stimulus.im_cell{ii}(:,:,1));
    title(stimulus.categories{ii});
end

% Compute power spectra for all categories ('im_cell')
D = loadDisplayParams('HiResDefault');

figure('Name', 'NEW STIM SF');hold on
peaks = [];
for ii = 1:24%length(stimulus.im_cell)
    
    im = stimulus.im_cell{ii};
    [frequencies, amplitudes, binnedFrequencies, binnedAmplitudes] = ...
        compute2DamplitudeSpectrum(im, D);
    [y,x] = max(binnedAmplitudes);
    peaks(ii) = binnedFrequencies(x);
    
    % Plot spectrum per category
    subplot(4,6,ii);
    plot(binnedFrequencies, binnedAmplitudes/max(binnedAmplitudes), 'LineWidth', 2);
    xlabel('Cycles per degree')
    xlim([0 7.5])
    set(gca, 'XGrid', 'on', 'XTick', 0:1.5:7.5)
    title(stimulus.categories{ii});
end

% Plot distribution of peaks
figure('Name', 'NEW STIM SF PEAKS'); 
subplot(1,2,1);histogram(peaks); xlabel('max spatial frequency'); ylabel('number of images')
xlim([0 4]);
subplot(1,2,2);plot(peaks, 'o-'); ylabel('category index'); ylabel('max spatial frequency'); 
xlim([0 25]);

% Compute luminance histograms

figure ('Name', 'NEW STIM LUM');hold on
for ii = 1:24%length(stimulus.im_cell)
    
    im = stimulus.im_cell{ii};
    % Plot histogram per category
    subplot(4,6,ii);
    histogram(im(im~=mode(im(:)))); title('(mode excluded from hist)');
    
    xlabel('Pixel values')
    xlim([0 255])
    title(stimulus.categories{ii});
end

%% DO THE SAME FOR THE ORIGINAL KAY STIMULI

% Pattern stimuli from Kay et al 2013
load('/Volumes/server/Projects/BAIR/Stimuli/Kay2013_stimuli/stimuli.mat');

% From http://kendrickkay.net/socmodel/index.html#contentsofstimuli

% Contents of 'stimuli.mat':
% 
% 'images' contains the raw stimulus frames for stimulus sets 1, 2, and 3
% (concatenated). The 'images' variable is a cell vector of dimensions 1 x
% 260 (since 69+156+35=260). Each entry corresponds to one stimulus, and
% each stimulus consists of one or more frames. For stimulus set 1 (1
% through 69), all entries consist of 30 distinct frames. For stimulus set
% 2 (70 through 225), all entries consist of 9 distinct frames except for
% entry #174 which consists of 7 distinct frames. For stimulus set 3 (226
% through 260), all entries consist of a single frame. The resolution of
% the images in stimulus set 1 is 600 pixels x 600 pixels; the resolution
% of the images in stimulus sets 2?3 is 800 pixels x 800 pixels. For all
% images, the format is uint8; the range of values is [0,254]; and the
% background has a value of 127.
% 
% 'conimages' contains the spatial masks (contrast images) used in the
% generation of some of the stimuli. There is a direct correspondence
% between 'conimages' and 'images'. The 'conimages' variable is a cell
% vector of dimensions 1 x 260. Each entry gives the mask that was used to
% generate the corresponding stimulus in 'images'. The indices of the
% stimuli that have associated masks are 1:69, 70:138, and 185:208. For the
% other indices, the entry in 'conimages' is simply the empty matrix. The
% resolution of the contrast images in stimulus set 1 is 600 pixels x 600
% pixels; the resolution of the contrast images in stimulus sets 2?3 is 256
% pixels x 256 pixels. For all contrast images, the format is double; and
% the range of values is [0,1] where a value of X indicates that at that
% pixel, the underlying stimulus, weighted by X, was blended with the gray
% background, weighted by 1?X.
% 
% 'bpfilter' contains the band-pass filter used to generate some of the
% stimuli in stimulus sets 2?3. The filter is a matrix of dimensions 21 x
% 21 and was used in stimulus construction, at which point the image
% resolution was 256 pixels x 256 pixels.


knk_idx = [...
            162 ... CRF-1
            164 ... CRF-2
            166 ... CRF-3
            167 ... CRF-4
            116 ... CRF-5
            150 ... GRATING
            154 ... PLAID
            158 ... CIRCULAR
            184 ... SPARSITY-1
            183 ... SPARSITY-2
            182 ... SPARSITY-3
            181 ... SPARSITY-4
            171 ... FACES-1
            171 ... FACES-2
            171 ... FACES-3
            171 ... FACES-4
            173 ... LETTERS-1
            173 ... LETTERS-2
            173 ... LETTERS-3
            173 ... LETTERS-4
            172 ... SCENES-1 (175??)
            172 ... SCENES-2 (175??)
            172 ... SCENES-3 (175??)
            172 ... SCENES-4 (175??)
            116 ... ONEPULSE-1
            116 ... ONEPULSE-2
            116 ... ONEPULSE-3
            116 ... ONEPULSE-4
            116 ... ONEPULSE-5
            116 ... ONEPULSE-6
            116 ... TWOPULSE-1
            116 ... TWOPULSE-2
            116 ... TWOPULSE-3
            116 ... TWOPULSE-4
            116 ... TWOPULSE-5
            116 ... TWOPULSE-6
            ];

whichIm = {...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:2 ...
            3:4 ...
            5:6 ...
            7:8 ...
            1:2 ...
            3:4 ...
            5:6 ...
            7:8 ...
            1:2 ...
            3:4 ...
            5:6 ...
            7:8 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            1:6 ...
            };

knk_im_cell = [];
for ii = 1:24%length(knk_idx)
    knk_im_cell{ii} = images{knk_idx(ii)}(:,:,whichIm{ii}(1));
end

figure('Name', 'OLD STIM');hold on
for ii = 1:24%length(knk_im_cell)
    subplot(4,6,ii);
    imshow(knk_im_cell{ii}(:,:,1));
    title(stimulus.categories{ii});
end

% D = loadDisplayParams('cni_lcd');
D = loadDisplayParams('cni_lcd_mock');

figure('Name', 'OLD STIM SF');hold on
peaks = [];
for ii = 1:24%length(knk_im_cell)
    
    im = knk_im_cell{ii}(:,:,1);
    [frequencies, amplitudes, binnedFrequencies, binnedAmplitudes] = ...
        compute2DamplitudeSpectrum(im, D);
    [y,x] = max(binnedAmplitudes);
    peaks(ii) = binnedFrequencies(x);
    
    % Plot spectrum per category
    subplot(4,6,ii);
    plot(binnedFrequencies, binnedAmplitudes/max(binnedAmplitudes), 'LineWidth', 2);
    xlabel('Cycles per degree')
    xlim([0 7.5])
    set(gca, 'XGrid', 'on', 'XTick', 0:1.5:7.5)
    title(stimulus.categories{ii});
end

% Plot distribution of peaks
figure('Name', 'OLD STIM SF PEAKS'); 
subplot(1,2,1);histogram(peaks); xlabel('max spatial frequency'); ylabel('number of images')
xlim([0 4]);
subplot(1,2,2);plot(peaks, 'o-'); ylabel('category index'); ylabel('max spatial frequency'); 
xlim([0 25]);

% Compute luminance histograms
figure('Name', 'OLD STIM LUM');hold on
for ii = 1:24%length(knk_im_cell)
    
    im = knk_im_cell{ii};
    % Plot histogram per category
    subplot(4,6,ii);
    histogram(im(im~=mode(im(:)))); title('(mode excluded from hist)');
    
    xlabel('Pixel values')
    xlim([0 255])
    title(stimulus.categories{ii});
end

%% make stimulus movie

% movie ----
%     movieName = sprintf('~/Desktop/spatiotemporal%02d.avi', runnum);
%     nFramePerSec = 60;
%
%     v = VideoWriter(movieName);
%     v.FrameRate = nFramePerSec;
%     v.Quality   = 100;
%     open(v)
%     c = 1;
%     %fH = figure(); set(fH, 'Visible', 'off')
%
%     cmap = gray(256);
%
%     fprintf('Making scan number %d\n', runnum);
%     for ii = 1 : length(stimulus.seq)
%
%
%         im = uint8(stimulus.images(:,:,stimulus.seq(ii)));
%         frame = im2frame(im, cmap);
%
%         writeVideo(v, frame);
%
%         if mod(ii,100) == 0, fprintf('.'); drawnow(); end
%
%     end
%     close(v)