%% paths
cd('/Volumes/server/Projects/BAIR/Stimuli/Code')
tbUse bairStimuli;
addpath(genpath('/Volumes/server/Projects/BAIR/Stimuli/Code/NYU ECOG WINAWERLAB CODE/vistadisp'));
 
% read stim params directly from
d = loadDisplayParams('SoMMacBook');

%% Spatiotemporal
load('/Volumes/server/Projects/BAIR/ECoG/648/stimdata/648_20171110T142403.mat')
for ii = 1:36
    figure, imshow(stimulus.im_cell{ii}(:,:,1));
    title(stimulus.categories{ii});
    set(gcf, 'numbertitle', 'off', 'Name', stimulus.categories{ii});
end

%% HRF
load('/Volumes/server/Projects/BAIR/ECoG/648/stimdata/648_20171110T142154.mat')
figure, imshow(stimulus.images(:,:,1))
set(gcf, 'numbertitle', 'off', 'Name', 'HRF');

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

%% compare with rotational average?
[Zr, R] = radialavg(A,length(binCenters));
figure;plot(R,Zr);
figure;plot(binCenters,Zr);

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
for ii = 1:size(stimulus.im_cell,2)
    
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
