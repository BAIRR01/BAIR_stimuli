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

%         % Download spatiotemporal stimuli 
%         images = [];
%         for ii = 1:13 % why 13?
%             readPth = sprintf('https://wikis.nyu.edu/download/attachments/85394548/spatiotemporal%d.mat?api=v2', ii);
%             stimDir = fullfile(BAIRRootPath, 'stimuli');
%             fname = sprintf('spatiotemporal%d.mat', ii);
%             writePth = fullfile(stimDir, fname);
%             if ~exist(writePth, 'file'),  websave(writePth,readPth); end
%             tmp = load(writePth);
%             images = [images tmp.im];
%         end 


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

        
        im_cell = images(knk_idx);
        
        images = [];
        for ii = 1:num_cats
            these_images = im_cell{ii}(:,:,whichIm{ii});
            images = cat(3, images, these_images);
        end
        
        
%% durations

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

%% indices into original unfitlered stim:
% LETTERS
%                             imageIndex = [13 21; 
%                                           6 20; 
%                                           1 9; 
%                                           8 24];
% HOUSES
%                             imageIndex = [33 7;  
%                                           24 20; 
%                                           41 3; 
%                                           31 25]; 
% FACES
%                             imageIndex = [30 35;
%                                           24 83; 
%                                           74 13; 
%                                           59 41];     

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
