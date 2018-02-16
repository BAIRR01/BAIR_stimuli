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

%% QUESTIONS

% GRATING / CIRCULAR / PLAID: % how to equate spacing of dark and light
%   lines? PLAID looks like pattern stimulus: maybe bp filter should be
%   adapted for these stimuli?

% FACES / HOUSES / LETTERS: resolution sufficient? do we want the exact
%   same images as in the SOC paper? are some images repeated?

%% Make the images

% determine if we're creating the master or loading & resizing for a specific display
site = stimParams.experimentSpecs.Row{1};
imageSizeInPixels = size(stimParams.stimulus.images);
      
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
        
        % Create CRF
        
        numberOfCategories = length(find(contains(categories, 'CRF')));
        numberOfImagesPerCat = 6;
        CRFimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        CRFim_cell = cell([numberOfCategories 1]);
        
       % From Kay et al, 2013 PLOS CB
%         These stimuli were constructed by varying the contrast of the
%         noise patterns used in SPACE. Ten different contrast levels were
%         used: 1%, 2%, 3%, 4%, 6%, 9%, 14%, 21%, 32%, and 50%. These
%         contrast levels are relative to the contrast of the patterns used
%         in SPACE, which is taken to be 100%.

        contrastLevels = [0.0625 0.125 0.25 0.5 1]; % or should we use log/linear?
        densityLevels = 20; % middle density stimulus for SPARSITY

        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createPatternStimulus(stimParams, densityLevels, contrastLevels(cc));
                
                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                CRFimages(:,:,imCount) = image8Bit;
                
            end
            CRFim_cell{cc} = CRFimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

        % Create SPARSITY 
        
        numberOfCategories = length(find(contains(categories, 'SPARSITY')));
        numberOfImagesPerCat = 6;
        SPARSITYimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        SPARSITYim_cell = cell([numberOfCategories 1]);
               
       % From Kay et al, 2013 PLOS CB
%         We generated noise patterns using cutoff frequencies of 2.8, 1.6,
%         0.9, 0.5, and 0.3 cycles per degree, and numbered these from 1
%         (smallest separation) to 5 (largest separation). The noise
%         patterns used in SPACE correspond to separation 4; thus, we only
%         constructed stimuli for the remaining separations 1, 2, 3, and 5.
%         The noise patterns occupied the full stimulus extent (no aperture
%         masking).
        
        contrastLevels = 1; % or should we use log/linear?
        densityLevels = [10 15 25 30]; % middle density is CRF full contrast stim

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
        
        % Create TEMPORAL
        
        % First and second for double pulse should be same stimulus (TO DO
        % check how this is determined in stim sequencing)
        
        numberOfCategories = length(find(contains(categories, 'PULSE')));
        numberOfImagesPerCat = 6;
        PULSEimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        PULSEim_cell = cell([numberOfCategories 1]);
        
        contrastLevels = 1; % max contrast
        densityLevels = 20; % middle density
        
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
        
        % Create GRATING, PLAID, CIRCULAR
        
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
        
        numberOfCategories = length(find(contains(categories, {'GRATING'})));
        numberOfImagesPerCat = 6;
        GRATINGimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        GRATINGim_cell = cell([numberOfCategories 1]);

        contrastLevels = 1; % max contrast
        cpd = 3;            % degrees
        orientation = pi/2; % horizontal
        
        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createGratingStimulus(stimParams, cpd, contrastLevels, (1:16)/pi);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                GRATINGimages(:,:,imCount) = image8Bit;

            end
            GRATINGim_cell{cc} = GRATINGimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

        % DEBUG  % compare with old version
        % load('~/Box Sync/BAIR/Stimuli/spatiotemporal_NYU-ECOG_101.mat')
        % figure;imshow(stimulus.images(:,:,31)); title('OLD');
        % figure;imshow(image8Bit); title('NEW');

        numberOfCategories = length(find(contains(categories, {'PLAID'})));
        numberOfImagesPerCat = 6;
        PLAIDimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        PLAIDim_cell = cell([numberOfCategories 1]);

        contrastLevels = 1; % max contrast
        densityLevels = 50; % middle density
        orientation = pi/4; 
        
        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createPlaidStimulus(stimParams, densityLevels, contrastLevels, orientation);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                PLAIDimages(:,:,imCount) = image8Bit;

            end
            PLAIDim_cell{cc} = PLAIDimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

%         DEBUG  % compare with old version
%         load('/Users/winawerlab/Box Sync/Stimuli/spatiotemporal_NYU-ECOG_101.mat')
%         figure;imshow(stimulus.images(:,:,37)); title('OLD');
%         figure;imshow(image8Bit); title('NEW');

        numberOfCategories = length(find(contains(categories, {'CIRCULAR'})));
        numberOfImagesPerCat = 6;
        CIRCULARimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        CIRCULARim_cell = cell([numberOfCategories 1]);

        contrastLevels = 1; % max contrast
        densityLevels = 10; % middle density
        
        for cc = 1:numberOfCategories
            for ii = 1:numberOfImagesPerCat
                imageForThisTrial = createCircularStimulus(stimParams, densityLevels, contrastLevels);

                % Double to unsigned 8 bit integer, needed for vistadisp
                image8Bit = uint8((imageForThisTrial+.5)*255);
                imCount = ii+(cc-1)*numberOfImagesPerCat;
                CIRCULARimages(:,:,imCount) = image8Bit;

            end
            CIRCULARim_cell{cc} = CIRCULARimages(:,:,(cc-1)*numberOfImagesPerCat+1);
        end

%         DEBUG  % compare with old version
%         load('/Users/winawerlab/Box Sync/Stimuli/spatiotemporal_NYU-ECOG_101.mat')
%         figure;imshow(stimulus.images(:,:,43)); title('OLD');
%         figure;imshow(image8Bit); title('NEW');

        % Create FACES / HOUSES / LETTERS
        
        numberOfCategories = length(find(contains(categories, {'FACES', 'LETTERS', 'HOUSES'})));
        numberOfImagesPerCat = 2;
        FLHimages = uint8(zeros([imageSizeInPixels numberOfCategories * numberOfImagesPerCat]));
        FLHim_cell = cell([numberOfCategories 1]);
     
%         % load orig stim from Kay et al., 2013
        load('/Volumes/server/Projects/BAIR/Stimuli/Kay2013_unfilteredstimuli/stim.mat');
         
        % - what you have for 'faces' and 'houses' are values between
        % [0,1].  They are intended to be displayed on standard monitors
        % (with squaring).  You may want to apply .^2 to the values to
        % simulate what luminance values are.
        % 
        % - The SOC stimuli involve only 9 instances of these classes.  I
        % didn't keep a record of the random number seed, so you'll have to
        % figure out which ones they were.
        % 
        % - There is some posthoc imresizing that I do to create the
        % official SOC stimuli...
        % 
        % - There may also be some masking, too (i.e. blending into the
        % gray background)

        % From Kay et al, 2013 PLOS CB
%          Each image was whitened (to remove low-frequency bias) and then
%          filtered with the custom band-pass filter (described
%          previously). Finally, the contrast of each image was scaled to
%          fill the full luminance range.
        
%         I = faces(:,:,1);
%          
%         % resize to match bp filter resolution
%         I1 = imresize(I, imageSizeInPixels); 
%         
%         % whiten (?)
%         % whiten = @(x) (x - mean(x(:)))./ diff(prctile(x(:), [.1 .9]));
%         I2 = I1 - mean(I1(:));
%         
%         % take square to 'simulate luminance' values (see above)
%         I3 = I2.^2;
%        
%         % band-pass filter
%         I4 = imfilter(I3, stimParams.bpFilter);
%         I5 = uint8((I4+.5)*255);
%         
%         figure;
%         subplot(2,3,1);imshow(I0);axis on; title('original');
%         subplot(2,3,2);imshow(I1);axis on; title('resized');
%         subplot(2,3,3);imshow(I2);axis on; title('whiten');
%         subplot(2,3,4);imshow(I3);title('squared');
%         subplot(2,3,5);imshow(I4);title('filter');
%         subplot(2,3,6);imshow(I5);title('8bit');
%         
        % why is the image (background) noisy? is it due to the resizing?
        % how about the contrast enhancement step in the SOC description?
        % NOTE: old face stim seem to be cropped + pasted into background
        
        % Concatenate all category types
        images = cat(3, CRFimages, GRATINGimages, PLAIDimages, CIRCULARimages, SPARSITYimages, FLHimages, PULSEimages);
        im_cell = cat(1, CRFim_cell, GRATINGim_cell, PLAIDim_cell, CIRCULARim_cell, SPARSITYim_cell, FLHim_cell, PULSEim_cell);
        
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
        
    otherwise    
        
        % Resize the Master stimuli to the required stimulus size for this
                % modality and display
        fprintf('[%s]: Loading Master stimuli for: %s\n', mfilename, site);

        % Load the Master stimuli
        master_stimulus = loadBAIRStimulus('spatiotemporal', 'Master', runNum);
        
        % Resize         
        fprintf('[%s]: Resizing Master stimuli for: %s\n', mfilename,  site);
        imageSizeInPixels = size(stimParams.stimulus.images);
        images = imresize(master_stimulus.images, imageSizeInPixels );
        
        im_cell = cell(num_cats);
        for ii = 1:num_cats
            im_cell{ii} = imresize(master_stimulus.im_cell{ii}, imageSizeInPixels);
        end
        
        % Soft circular mask (1 pixel of blurring per 250 pixels in the image)
        supportDiameter       = imageSizeInPixels(1);
        maskRadius            = (stimParams.stimulus.srcRect(3) - stimParams.stimulus.srcRect(1))/2;
        circularMask          = mkDisc(supportDiameter, maskRadius, (imageSizeInPixels+1)./2, 1/250 * imageSizeInPixels(1));
        imagesDouble          = double(images)/255-.5;
        imagesMasked          = bsxfun(@times,imagesDouble,  circularMask);
        images                = uint8((imagesMasked+.5)*255);
        
        % copy these from master 
        categories = master_stimulus.categories;
        durations = master_stimulus.duration;
        ISI = master_stimulus.ISI;
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

num_cats = length(categories);
ITIs = linspace(ITI_min,ITI_max,num_cats);

% This is the stimulus structure used by vistadisp
stimulus.cmap         = stimParams.stimulus.cmap;
stimulus.srcRect      = stimParams.stimulus.srcRect;
stimulus.dstRect      = stimParams.stimulus.destRect;

stimulus.im_cell      = im_cell;
stimulus.images       = images;

stimulus.categories   = categories;
stimulus.duration     = durations;

stimulus.ISI          = ISI;
stimulus.ITI          = ITIs;
stimulus.prescan      = prescan; % seconds
stimulus.postscan     = postscan; % seconds

% make individual trial sequences
%   randomize stimulus order

% TO DO: update according to new stim generation
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
diodeSeq = zeros(size(stimulus.seq));

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

fname = sprintf('spatiotemporal_%s_%d', site, runNum);

save(fullfile(vistadispRootPath, 'Retinotopy', 'storedImagesMatrices',  fname), 'stimulus', '-v7.3')

return
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

%% DEBUG : plot individual stim
whichCat = 'PULSE';

images = eval([whichCat 'images']);
im_cell = eval([whichCat 'im_cell']);
category_index = find(contains(categories, whichCat));

figure;hold on
for ii = 1:length(im_cell)
    subplot(1,length(category_index),ii);
    imshow(im_cell{ii}(:,:,1));
    title(categories{category_index(ii)});
end

figure;hold on
for ii = 1:size(images,3)
    subplot(length(category_index),size(images,3)/length(category_index),ii);
    imshow(images(:,:,ii));
end

%% DEBUG : plot all exemplar stim, compute SF

%% Plot all stim
figure;hold on
for ii = 1:length(stimulus.im_cell)
    subplot(6,6,ii);
    imshow(stimulus.im_cell{ii}(:,:,1));
    title(stimulus.categories{ii});
end

% read stim params directly from
d = loadDisplayParams('HiResDefault');

%% Power spectrum

sz = size(stimulus.images(:,:,1));

% compute bins from first image
stimSz = 2*pix2angle(d,sz(1)/2);
stim = double(stimulus.images(:,:,1));
stim = stim-mode(stim(:));
A = abs(fftshift(fft2(stim))); % can pad by putting extra args
Fs1D = (-sz(1)/2:sz(1)/2-1)/stimSz;

[FsX, FsY] = meshgrid(Fs1D);
F = sqrt(FsX.^2+FsY.^2);

[~,edges,bin] = histcounts(F(:));
binCenters = (edges(1:end-1)+edges(2:end))/2;

figure, imagesc([min(FsX) max(FsX)] , [min(FsY) max(FsY)], abs(A))
Mn = accumarray(bin, A(:), [], @mean);
figure, plot(F(:), A(:), '.k', binCenters, Mn)

%% Loop across images to extract bin average for each image, plot 

figure, hold on

peaks = zeros(size(stimulus.images,3),1);
for ii = 1:size(stimulus.images,3)
    
    stim = double(stimulus.images(:,:,ii));
    stim = stim-mode(stim(:)); 
    A = fftshift(abs(fft2(stim)));

    Mn = accumarray(bin, A(:), [], @mean);
    plot(binCenters, Mn)
    [y,x] = max(Mn);
    peaks(ii) = binCenters(x);
end

figure; hist(peaks); xlabel('max spatial frequency'); ylabel('number of images')


%% SOC by category

sz = size(stimulus.im_cell{1}(:,:,1));%ssize(stimulus.images(:,:,1));

% compute bins from first image
stimSz = 2*pix2angle(d,sz(1)/2);
stim = double(stimulus.im_cell{1}(:,:,1));%double(stimulus.images(:,:,1));
stim = stim-mode(stim(:));
A = fftshift(abs(fft2(stim)));
Fs1D = (-sz(1)/2:sz(1)/2-1)/stimSz;

[FsX, FsY] = meshgrid(Fs1D);
F = sqrt(FsX.^2+FsY.^2);

[N,edges,bin] = histcounts(F(:));
binCenters = (edges(1:end-1)+edges(2:end))/2;

figure, hold on
peaks = zeros(size(stimulus.im_cell,2),1);
for ii = 1:size(stimulus.im_cell,1)
    
    stim = double(stimulus.im_cell{ii}(:,:,1));
    stim = stim-mode(stim(:)); 
    A = fftshift(abs(fft2(stim)));

    Mn = accumarray(bin, A(:), [], @mean);
    subplot(6,6,ii)
    plot(binCenters, Mn, 'LineWidth', 2)
    title(stimulus.categories{ii});
    set(gca, 'XLim', [0 5], 'YLim', [0 10^5.5]);
    [y,x] = max(Mn);
    peaks(ii) = binCenters(x);
    %axis tight
end

