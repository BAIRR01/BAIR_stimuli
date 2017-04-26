function stimMakeSpatiotemporalExperiment(s_example, modality)
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

imsize = s_example.stimulus.srcRect(4);

% Download spatiotemporal stimuli
images = [];
for ii = 1:13
    readPth = sprintf('https://wikis.nyu.edu/download/attachments/85394548/spatiotemporal%d.mat?api=v2', ii);
    stimDir = fullfile(BAIRRootPath, 'stimuli');
    fname = sprintf('spatiotemporal%d.mat', ii);
    writePth = fullfile(stimDir, fname);
    if ~exist(writePth, 'file'),  websave(writePth,readPth); end
    tmp = load(writePth);
    images = [images tmp.im];
end

switch(lower(modality))
    case 'fmri'
        numruns = 4;
        ITI_min  = 3;
        ITI_max  = 6;
        prescan  = 9; % seconds
        postscan = 9; % seconds
    case {'ecog' 'eeg' 'meg'}
        numruns = 12;
        ITI_min  = 1.25;
        ITI_max  = 1.75;
        prescan  = 2; % seconds
        postscan = 2; % seconds
    otherwise
        error('Unknown modality')
end

for runnum = 1:numruns
    
    stimulus = [];
    
    
    num_cats = 36;
    
    switch lower(modality)
        case 'fmri'
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
    stimulus.dstRect = [128 0 896 768];
    %       images: [768×768×91 uint8]
    %          cmap: [256×3 double]
    %           seq: [182×1 double]
    %     seqtiming: [182×1 double]
    %        fixSeq: [182×1 double]
    %       srcRect: [0 0 768 768]
    %      destRect: [128 0 896 786]
    %       trigSeq: [1×182 double]
    %      diodeSeq: [182×1 double]
    
    
    
    fname = sprintf('spatiotemporal_%s_%d', modality, runnum);
    save(fullfile(BAIRRootPath, 'stimuli', fname), 'stimulus')
    
end


