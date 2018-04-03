%% Increase contrast in hRF stimulis

load('/Volumes/server/Projects/BAIR/Stimuli/Code/all stimuli/hrf_ECOG_3.mat');
stimnew = double(stimulus.images);
stimnew = stimnew - 128;
stimnew = stimnew *2;
stimnew = uint8(stimnew + 128);
figure; histogram(stimulus.images); title('hrf_ECOG1')


stimulus.images = stimnew;


save('/Volumes/server/Projects/BAIR/Stimuli/Code/NYU ECOG WINAWERLAB CODE/vistadisp/Retinotopy/storedImagesMatrices/hrf_ECOG_13.mat', 'stimulus');

load('/Volumes/server/Projects/BAIR/Stimuli/Code/NYU ECOG WINAWERLAB CODE/vistadisp/Retinotopy/storedImagesMatrices/hrf_ECOG_13.mat');
figure; histogram(stimulus.images); title('hrf_ECOG11')

%% Increase contrast in spatiotemporal stimuli
pth = '/Volumes/server/Projects/BAIR/Stimuli/Code/all stimuli/';
for jj = 1:12
    load(fullfile(pth, sprintf('spatiotemporal_ECOG_%d.mat', jj)));
    
    scale_factor = 1.5;
    clip_lim = [127 - 127/scale_factor, 127+127/scale_factor];
    
    
    enhancecontrast = @(im, scale_factor) uint8((double(im)-127)*scale_factor+127);
    
    scale_it(1:168)=1.5;
    scale_it(31:48)=2;
    stimnew = stimulus.images;
    for ii = 1:168
        stimnew(:,:,ii) = enhancecontrast(stimulus.images(:,:,ii), scale_it(ii));
    end
    stimulus.images = stimnew;
    
    
    scale_it = ones(1,36)*1.5;
    scale_it(6:8)=2;
    stim_new = stimulus.im_cell;
    for ii = 1:36
        stim_new{ii} = enhancecontrast(stimulus.im_cell{ii}, scale_it(ii));
    end
    stimulus.im_cell = stim_new;
    
    %save(fullfile(pth, sprintf('spatiotemporal_ECOG_%d.mat', jj+100)), 'stimulus')

end
%%

figure,
for ii = 1:36
    subplot(4,9,ii)
    imshow((stimulus.im_cell{ii}(:,:,1)))
end

figure,
for ii = 1:36
    subplot(4,9,ii)
    idx = stimulus.im_cell{ii} ~= 127;
    %histogram(enhancecontrast(stimulus.im_cell{ii}(idx)))
    histogram(stimulus.im_cell{ii}(idx))

    hold on;
    %plot(clip_lim(1)*[1 1], get(gca, 'YLim'), 'k--');
    %plot(clip_lim(2)*[1 1], get(gca, 'YLim'), 'k--');
    set(gca, 'XLim', [0 255])
end