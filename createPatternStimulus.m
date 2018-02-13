function [output, edge, thresh, result] = createPatternStimulus(stimParams, relCutoff, scaleContrast)
% CREATE PATTERN STIMULUS
%   sz - the desired image size
%   relCutoff - Define a relative cutoff in terms of the available frequencies
%   bpfilter - the convolutional bandpass filter to use, in space domain

if nargin < 3 || isempty(scaleContrast)
    % Scale the edge contrast (this value was found empirically to maximize
    % contrast whilst preventing too much clipping for SOC stimuli)
    scaleContrast = 0.18;
else
    scaleContrast = scaleContrast * 0.18;
end

stimWidth  = stimParams.stimulus.srcRect(3)-stimParams.stimulus.srcRect(1);
stimHeight = stimParams.stimulus.srcRect(4)-stimParams.stimulus.srcRect(2);

% Create a random seed
im = randn(stimHeight, stimWidth);
im = im./(max(im(:) - min(im(:)))) + 0.5;

% Find the midpoint pixels
mid = ceil((size(im)+1)/2);

% Create a soft round filter in Fourier space
radius = relCutoff;%  size(im,1)/2;
mask = mkDisc(size(im), radius, mid, radius/5);

% Filter the image
dft = fftshift(fft2(im));
mdft = mask.*dft;
result = ifft2(ifftshift(mdft));

% Threshold the filtered noise
thresh = result - min(result(:)) > (max(result(:)) - min(result(:)))/2;

% Grab edges with derivative filter
edge1 = [0, 0, 0; 1, 0, -1; 0, 0, 0];
edge2 = [0, 1, 0; 0, 0, 0; 0, -1, 0];
edge = -1*(imfilter(double(thresh), edge1, 'circular').^2 + imfilter(double(thresh), edge2, 'circular').^2);

% Scale the edge contrast
edge = scaleContrast*edge; 

% Filter convolutionally with bpfilter in the image domain
output = imfilter(edge, stimParams.bpFilter, 'circular');

return

end

%% Debug 

function checkKayStimuli
% Pattern stimuli from Kay et al 2013
load('/Users/jonathanwinawer/Downloads/stimuli.mat')
D = loadDisplayParams('cni_lcd');
%
figure; % set(gca, 'ColorOrder', parula(4)); hold on

stims = [181:183 85 184];% 176:180 ];
counter = 0;
for ii = stims
    counter = counter+1;
    subplot(2, length(stims),counter)
    im =  images{ii}(:,:,1);
    imshow(im);    title(ii)

end


for ii = stims
    subplot(2, length(stims),counter+length(stims))

    [frequencies, amplitudes, binnedFrequencies, binnedAmplitudes] = ...
        compute2DamplitudeSpectrum(im, D);
    %plot(frequencies(:), amplitudes(:), '-', binnedFrequencies, binnedAmplitudes);
    plot(binnedFrequencies, binnedAmplitudes/max(binnedAmplitudes));
    xlabel('Frequency (cycles per degree)')
    xlim([0 10])
    set(gca, 'XGrid', 'on', 'XTick', 0:1.5:50)
end

%legend(stimulus.categories(stims));
%xlim([0 6])
%%
figure
counter = 0;
for ii = 10:-2:1
    counter = counter + 1;
    subplot(1,5,counter)
    relCutoff = 1/(ii*5);
    im = createPatternStimulus(stimParams, relCutoff);
    imshow(im+.5)
    title(round(1/relCutoff));
    axis([0 size(im,1) 0 size(im,1)] *12.5/16.6)
end

end
