function [frequencies, amplitudes, binnedFrequencies, binnedAmplitudes] = ...
    compute2DamplitudeSpectrum(im, displayParams)

% enforce double precision
%im = double(im);
im = (double(im)./255)-0.5;

sz = size(im);
stimSz = 2*pix2angle(displayParams,sz(1)/2);
Fs1D = (-sz(1)/2:sz(1)/2-1)/stimSz;
[FsX, FsY] = meshgrid(Fs1D);
frequencies = sqrt(FsX.^2+FsY.^2);

% subtract the mean
im = im - mean(im(:));

% Use a window?
% im = im .* window2(sz(1), sz(2),'hann');
% im = im .* window2(sz(1), sz(2), @hann);

% Fourier transform
amplitudes  = fftshift(abs(fft2(im)));

[~,edges,bin] = histcounts(frequencies(:));
binnedFrequencies = (edges(1:end-1)+edges(2:end))/2;

binnedAmplitudes = accumarray(bin, amplitudes(:), [], @mean);

end