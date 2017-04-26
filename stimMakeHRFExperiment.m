function stimMakeHRFExperiment(s_example, modality)
%stimMakeHRFExperiment(s_example, modality)
%  Stimuli are presented for 0.25 seconds, with an exponentially
%  distributed interstimulus interval (mean ~9s, range 3-24s).

imsize = s_example.stimulus.srcRect(4);

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


lags = round(4*y);

hgexport(gcf, fullfile(BAIRRootPath, 'figures', 'HRF_ISIs.eps'))

seq = randperm(n);
onsets = cumsum([1 lags(seq)]);
figure(2), clf; set(gcf, 'Color', 'w')
set(gca, 'FontSize', 24, 'XTick', 0:60:300, 'YTick', []); hold on
stem(onsets, ones(1,n+1), 'LineWidth', 2)
xlabel('Time (s)')
hgexport(gcf, fullfile(BAIRRootPath, 'figures', 'HRF_onsets.eps'))


bar_carrier = stimMakeBarCarrier();

numruns = 12;

numim   = n+1;
images  = zeros(imsize, imsize, numim, 'uint8')+127;

bar_carrier = bar_carrier-0.5;
[x,y] = meshgrid(linspace(-1,1,imsize));
[th,r] = cart2pol(x,y);
mask = r<=1;


for ii = 1:numim-1
    images(:,:,ii) =   (bar_carrier(:,:,ii) .* double(mask))*127+128;
end




for runnum = 1:numruns
    seq    = randperm(n);
    onsets = cumsum(lags(randperm(n)));
    
    
    stimulus = [];
    stimulus.cmap       = s_example.stimulus.cmap;
    stimulus.srcRect    = s_example.stimulus.srcRect;
    stimulus.dstRect    = s_example.stimulus.dstRect;
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
    
    
    
    fname = sprintf('hrf_%s_%d', modality, runnum);
    save(fullfile(BAIRRootPath, 'Stimuli', fname), 'stimulus')
    subplot(4,4,runnum),
    plot(stimulus.seqtiming, stimulus.seq)
end
