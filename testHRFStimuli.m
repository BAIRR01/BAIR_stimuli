% % Fom Kay et al, 2013:
% Stimuli consisted of grayscale images restricted to a band-pass
% range of spatial frequencies centered at 3 cycles per degree. To
% enforce this restriction, a custom band-pass filter was used in the
% generation of some of the stimuli. The filter was a zero-mean
% isotropic 2D Difference-of-Gaussians filter whose amplitude
% spectrum peaks at 3 cycles per degree and drops to half-maximum
% at 1.4 and 4.7 cycles per degree. Restricting the spatial frequency
% content of the stimuli avoids the complications of building multiscale
% models and helps constrain the scope of the modeling
% endeavor. Even with the spatial frequency restriction, it is possible
% to construct a rich diversity of stimuli including objects and other
% naturalistic stimuli.

% Get stimParams from s_Make_BAIR_Visual_Experiments

% Load a stored bandpass filter
load('/Users/jonathanwinawer/matlab/junk/bpfilter.mat')
bpfilter = imresize(bpfilter, 2); % 42 x 42 pixels


% pixels in stimulus diameter: 728
stimDiameterInPixels = stimParams.stimulus.srcRect(4);
% % Make a bandpass filter from a laplacian of gaussian
% %bpfilter = -fspecial('log', 50, 4);
% 
% Make a bandpass filter from a difference of two gaussians
filterSz       = round(stimDiameterInPixels/42);
filterCenter   = fspecial('gaussian', filterSz, filterSz/15.5); % high frequency cut-off
filterSurround = fspecial('gaussian', filterSz, filterSz/10.5);  % low  frequency cut-off
%bpfilter =  filterCenter/(mean(filterCenter(:))) - filterSurround/mean(filterSurround(:));


bpfilter = bpfilter / norm(bpfilter(:));

subplot(1,2,1)
surf(bpfilter)

sz = size(bpfilter);
stimSz = 2*pix2angle(stimParams.display,sz(1)/2);
Fs1D = (-sz(1)/2:sz(1)/2-1)/stimSz;
[FsX, FsY] = meshgrid(Fs1D);
F = sqrt(FsX.^2+FsY.^2);

A  = fftshift(abs(fft2(bpfilter)));


[~,edges,bin] = histcounts(F(:));
binCenters = (edges(1:end-1)+edges(2:end))/2;

Mn = accumarray(bin, A(:), [], @mean);
subplot(1,2,2); hold on
plot(F(:), A(:), '-');%, binCenters, Mn);
set(gca, 'XTick', [1.4 3 4.7], 'XGrid', 'on', 'XLim', [0 10])
xlabel('Frequency (cycles per degree)')


