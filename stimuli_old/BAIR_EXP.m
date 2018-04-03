%% QUESTIONS FOR UTRECHT
% HOW DO YOU DEAL WITH TRIGGERS (FMRI)/ ECOG?
% CALIBRATION FILES (INC SCREEN RESOLUTION)
% HOW DO YOU RECEIVE KEY PRESSES? (65-68 for fMRI)
% WHAT SOFTWARE DO YOU NORMALLY USE FOR PRFS?

imsize = 768;
% %% CHECK ON GAMMA TABLES INC NYU SOM, UTRECHT IEMU, 3T/7T
% 
% t = 0:.01:216;
% s = zeros(size(t));
% s(mod(t,24)<12) = 1;
% figure, set(gcf, 'Color', 'w')
% plot(t,s)
% axis([-5 221 -.1 1.1])
% set(gca, 'YTick', [0 1], 'YTickLabel', {'Rest', 'Task'},'XTick', 0:24:216, 'FontSize', 18)
% hgexport(gcf, '~/Desktop/task.eps')

%% MAKE TASK EXPERIMENT
modality = 'MEG';
a=load(sprintf('spatiotemporal_%s_1.mat', modality));

switch(modality)
    case 'fMRI'
        numruns = 2;
    case 'MEG'
        numruns = 2;
end

for runnum = 1:numruns
    stimulus = [];
    stimulus.cmap       = a.stimulus.cmap;
    stimulus.srcRect    = a.stimulus.srcRect;
    stimulus.dstRect    = a.stimulus.destRect;
    stimulus.images     = a.stimulus.images(:,:,end);
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
    
    pth = '~/matlab/toolboxes/BAIR/stimuli/';
    fname = sprintf('task_%s_%d', modality, runnum);
    save(fullfile(pth, fname), 'stimulus')
    
end
% Anticipatory BOLD
% runnum = 1;
% params = retCreateDefaultGUIParams;
% params.fixation      = 'contrastdecrementdot';
% params.skipSyncTests = false;
% params.triggerKey    = '`'; % comes with each volume (confirm which key for Utrecht)
% params.calibration   = 'CBI_NYU_projector';
% params.experiment    = 'experiment from file';
% loadMatrix_str       = sprintf('antBOLDnoStimulus%d.mat', runnum);
% params.loadMatrix    = loadMatrix_str;
% params.modality      = 'fMRI';
%
% ret(params)

%% HRF (zebras)
%  Stimuli are presented for 0.25 seconds, with an exponentially
%  distributed interstimulus interval (mean ~9s, range 3-24s).
n    = 32;
mint = 3;
maxt = 24;
x = linspace(0,1,n) * (1-exp(-mint)) + exp(-mint);
y= -log(x)/mint;
y = y*(maxt-mint)+mint;
disp(mean(y))
figure(1),clf, set(gcf, 'Color', 'w')
set(gca, 'FontSize', 24); hold on
plot(y, 'o-', 'LineWidth', 4, 'MarkerSize', 12); axis tight
ylabel('ITI (s)'); xlabel('Trial')
%hgexport(gcf, '~/Desktop/HRF_ISIs.eps')

lags = round(4*y);

a=load('/Users/winawerlab/matlab/git/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices/spatiotemporal_MEG_12.mat');
numruns = 12;
numim   = n+1;
images  = zeros(768, 768, numim, 'uint8')+127;

load('bar_carrier.mat')

bar_carrier = bar_carrier-0.5;
[x,y] = meshgrid(linspace(-1,1,imsize));
[th,r] = cart2pol(x,y);
mask = r<=1;


for ii = 1:numim-1;
    images(:,:,ii) =   (bar_carrier(:,:,ii) .* double(mask))*127+128;
end
figure, 
for runnum = 1:numruns
    seq    = randperm(n);
    onsets = cumsum(lags(randperm(n)));
    
    
    stimulus = [];
    stimulus.cmap       = a.stimulus.cmap;
    stimulus.srcRect    = a.stimulus.srcRect;
    stimulus.dstRect    = a.stimulus.destRect;
    stimulus.images     = images;
    stimulus.seqtiming  = (1:300*4)/4;
    stimulus.seq        = zeros(size(stimulus.seqtiming))+numim;
    stimulus.seq(onsets) = seq;
    stimulus.fixSeq     = ones(size(stimulus.seqtiming));    

    this_frame = 0;
    while true
        % wait between 4 and 20 frames (1 to 5 seconds)
        isi = randperm(16,1)+3;
        this_frame = this_frame + isi;
        if this_frame > length(stimulus.fixSeq), break; end
        stimulus.fixSeq(this_frame:end) = mod(stimulus.fixSeq(this_frame-1),2)+1;
    end
    
    switch lower(modality)
        case 'fmri'
        otherwise
            stimulus.trigSeq  = double(stimulus.seq>0);
            stimulus.diodeSeq = stimulus.trigSeq;
    end
    
    
    pth = '~/matlab/toolboxes/BAIR/stimuli/';
    fname = sprintf('hrf_%s_%d', modality, runnum);
    save(fullfile(pth, fname), 'stimulus')
    subplot(4,4,runnum),
    plot(stimulus.seqtiming, stimulus.seq)
end


%% PRF
modality = 'fMRI';
a=load(sprintf('spatiotemporal_%s_1.mat', modality));
numruns = 2;

load('bar_apertures.mat')
load('bar_carrier.mat')

bar_carrier = bar_carrier-0.5;

numim = size(bar_apertures,3)*3;
% 3 images per aperture
images = zeros(768, 768, numim, 'uint8')+127;
for ii = 1:numim
    idx = ceil(ii/3);
    idx2 = randsample(size(bar_carrier,3),1);
    images(:,:,ii) = bar_apertures(:,:,idx) .* bar_carrier(:,:,idx2) * 127+127;
end



for runnum = 1:numruns
    stimulus = [];
    stimulus.cmap       = a.stimulus.cmap;
    stimulus.srcRect    = a.stimulus.srcRect;
    stimulus.dstRect    = a.stimulus.destRect;
    stimulus.images     = images;
    stimulus.seqtiming  = (0:numim-1)/2;
    stimulus.seq        = 1:length(stimulus.seqtiming);
    stimulus.fixSeq     = ones(size(stimulus.seqtiming));    
    
   
    
    pth = '~/matlab/toolboxes/BAIR/stimuli/';
    fname = sprintf('ret_%s_%d', modality, runnum);
    save(fullfile(pth, fname), 'stimulus')
    
end

% HRF
%       Exponential ISI, mean XX, range XX
% TASK
%       12 s fixation task, 12 s relax, repeat 8 times (to measure anticipatory BOLD)
% RETINOTOPY
%       Like Dumoulin and Wandell, 2008: 8 sweeps ? 4 cardinal, 4 diagonal (diagonals include 50% blanks)
% SPATIOTEMPORAL
%        VISUAL: 36 unique stimuli, shown once each per scan (0.5 s except for temporal stimuli),
%                with mean ISI of 4.5 s, range 3-6 s; orientation (3; 1 grating, 1 plaid, 1 circular);
%                contrast (5; noise patterns); spacing: (5: noise patterns, 1 overlaps with contrast);
%                objects (12: 3 faces, 3 scenes, 3 objects, 3 bodies); temporal (12; 6 durations; 6 ISIs);


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


% make stimulus struct
load('/Volumes/server/Projects/SOC/knk/socmodel/stimuli.mat');


%% SPATIOTEMPORAL
modality = 'MEG';
for runnum = 2:12
    
    stimulus = [];
    
    
    num_cats = 36;
    
    switch lower(modality)
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
        'SCENES-1' ...
        'SCENES-2' ...
        'SCENES-3' ...
        'SCENES-4' ...
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
    
    tmp = cumsum(cellfun(@length, whichIm));
    whichIm_Idx = [[0 tmp(1:end-1)]+1; tmp];
    
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
    
    
    ITIs = linspace(ITI_min,ITI_max,num_cats);
    
    stimulus.im_cell     = images(knk_idx);
    stimulus.images = [];
    for ii = 1:num_cats
        these_images = stimulus.im_cell{ii}(:,:,whichIm{ii});
        stimulus.images = cat(3, stimulus.images, these_images);
    end
    stimulus.images(:,:,end+1) = mode(stimulus.images(:));
    BLANK = size(stimulus.images,3);
    
    % resize
    stimulus.images = imresize(stimulus.images, [imsize, imsize]);
    
    stimulus.categories = categories;
    stimulus.duration   = durations;
    stimulus.ISI        = ISI;
    stimulus.ITI        = ITIs;
    stimulus.prescan    = prescan; % seconds
    stimulus.postscan   = postscan; % seconds
    
    
    % make individual trial sequences
    %   randomize stimulus order
    stim_seq = randperm(num_cats);
    for ii = 1:num_cats
        idx = stim_seq(ii);
        
        % choose exemplar based on run number
        possible_exemplars = whichIm_Idx(1,idx):whichIm_Idx(2,idx);
        n = length(possible_exemplars);
        this_exemplar = mod(runnum, n);
        
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
    switch lower(modality)
        case 'fmri'
        otherwise
            stimulus.trigSeq = trigSeq;
            stimulus.diodeSeq = diodeSeq;
    end
    
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
    
    stimulus.cmap = (1:256)' * ones(1,3);
    
    stimulus.srcRect = [0 0 imsize imsize];
    stimulus.destRect = [128 0 896 768];
    %       images: [768×768×91 uint8]
    %          cmap: [256×3 double]
    %           seq: [182×1 double]
    %     seqtiming: [182×1 double]
    %        fixSeq: [182×1 double]
    %       srcRect: [0 0 768 768]
    %      destRect: [128 0 896 786]
    %       trigSeq: [1×182 double]
    %      diodeSeq: [182×1 double]
    
    
    pth = '~/matlab/toolboxes/BAIR/stimuli/';
    fname = sprintf('spatiotemporal_%s_%d', modality, runnum);
    if ~exist(pth, 'dir'), mkdir(pth); end
    save(fullfile(pth, fname), 'stimulus')
    
end



%%

% test it
figure
for ii = 1:length(stimulus.seq_sparse)-1
    imshow(stimulus.images(:,:,stimulus.seq_sparse(ii)));
    axis image; colormap gray;
    pause(stimulus.seqtiming_sparse(ii+1)-stimulus.seqtiming_sparse(ii));
end

stimulus.cmap = (ones(3,1)*(1:256))';
stimulus.srcRect = [0 0 size(stimulus.images,1) size(stimulus.images,1)];

if saveStimulus
    stimulus = rmfield(stimulus, 'im_cell');
    save ~/matlab/BAIR_spatiotemporal stimulus
end

% FIXATION
% MULTIPLE RUNS


%% Run it

params = retCreateDefaultGUIParams;
params.experiment = 'experiment from file';
params.loadMatrix = 'BAIR_spatiotemporal.mat';
params.calibration = 'gunjou';
ret(params)


%% HRF
x = 3:.01:24;

y = exp(-(x-3)/21);

x = 0:.001:1;
y = exp(-x);
% log(y) = -(x-3)/21;
% 21*log(y) = -(x-3);
% -21*log(y)+3 = x;

figure(1); clf;
plot(x,y)

sum(x.*y)

%%
x = rand(1,100000);
y = -log(x);
y(y>4) = 4;
y = y/4;
y = y*21+3;
figure(1), histogram(y)

%%
% log(y) = -(x-3)/21;
% 21*log(y) = -(x-3);
x = rand(1,100000) * (1-exp(-3)) + exp(-3);
y= -log(x)/3;
y = y*21+3;
disp(mean(y))
figure(1),
clf

histogram(y, 'Normalization', 'pdf')
hold on

x = 0:.001:1;
y = exp(-x*3);
x = x*21+3;
y = y / (mean(y)*21);
plot(x,y)