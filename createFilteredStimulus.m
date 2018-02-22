function output = createFilteredStimulus(stimParams, unfilteredImage, imageProcessingParams)

% CREATE FILTERED STIMULUS
%   stimParams - defined using stimInitialize.m; should contain the
%       convolutional bandpass filter to use, in space domain
%       (stimParams.bpfilter)
%   unfilteredImage - The image to be bandpass filtered
%   imageProcessingParams - a struct that determines various settings for
%       resizing, masking and contrast boosting the image (see
%       StimMakeSpatiotemporalExperiment.m for explanations)
%
%   Note: functions makes use of knkutils toolbox function
%    'makecircleimage' (for applying a square mask)
%   Note: we currently do not perform the prewhitening step described in
%    Kay 2013

imageSizeInPixels = size(stimParams.stimulus.images);

I0 = unfilteredImage; 

% Resize to a certain percentage of imageSizeInPixels
I1 = imresize(I0, imageProcessingParams.imageScaleFactor * imageSizeInPixels); 

stimSizeInPixels = size(I1);

% Paste into a background of size imageSizeInPixels
I2 = padarray(I1,[(imageSizeInPixels(1)-size(I1,1))/2 (imageSizeInPixels(1)-size(I1,1))/2], 0);

% Apply band-pass filter
I3 = imfilter(I2, stimParams.bpFilter);

% Apply mask to get rid of edges around stim and background noise
maskRadius             = imageProcessingParams.maskScaleFactor*stimSizeInPixels(1)/2;
squareMask             = makecircleimage(imageSizeInPixels(1),maskRadius,[],[], imageProcessingParams.maskSoftEdge*maskRadius,1) .* ...
                            makecircleimage(imageSizeInPixels(1),maskRadius,[],[], imageProcessingParams.maskSoftEdge*maskRadius,2);
I4                     = bsxfun(@times,I3,squareMask);

% TO DO: add possibility additionally apply a circular mask for the faces?
%maskOrigin             = [(imageSizeInPixels(1)+1)./2 (imageSizeInPixels(1)+1)./2];
%maskRadius             = 0.8*stimSizeInPixels(1)/2;
%circularMask           = makecircleimage(imageSizeInPixels(1),maskRadius,[],[],1.1*maskRadius,0, [maskOrigin(1)+0.1*maskOrigin(1) maskOrigin(2)], [1 1.2]);
%combinedMask           = circularMask .* squareMask;
%I4                     = bsxfun(@times,I3,combinedMask);

% Boost contrast, but prevent (too much) clipping
output = I4 * imageProcessingParams.contrastRange /range(I4(I4~=0));

return

%% DEBUG

% use a more robust range for contrast boosting?
% scaleFactor = range(prctile(I4(I4~=0), [2 98]));

% Plot the results of each processing step

% Conversion to 8Bit (this happens outside of this function in
% StimMakeSpatiotemporalExperiment, but we do want do it here for plotting)
I5 = output;
I6 = uint8((I5+.5)*255);

figure;
subplot(3,4,1);imshow(I0);axis on; title('original');
subplot(3,4,2);histogram(I0);
set(gca, 'XLim', [0 1], 'YLim', [0 3*10^4]);

subplot(3,4,3);imshow(I2);axis on; title('resize & fit in master display size');
subplot(3,4,4);histogram(I2);
set(gca, 'XLim', [0 1], 'YLim', [0 2*10^4]);

subplot(3,4,5);imshow(I3);axis on; title('bp filter');
subplot(3,4,6);histogram(I3);
set(gca, 'XLim', [-5 5], 'YLim', [0 2*10^4]);

subplot(3,4,7);imshow(I4);axis on; title('circular mask');
subplot(3,4,8);histogram(I4);
set(gca, 'XLim', [-5 5], 'YLim', [0 2*10^4]);

subplot(3,4,9);imshow(I5);axis on;title('scale (to prevent too much clipping)');
subplot(3,4,10);histogram(I5(I5~=mode(I5(:)))); title('(mode excluded from hist)');
set(gca, 'XLim', [-1 1])% 'YLim', [0 0.1*10^4]);

subplot(3,4,11);imshow(I6);axis on;title('8bit');
subplot(3,4,12);histogram(I6(I6~=mode(I6(:)))); title('(mode excluded from hist)');
set(gca, 'XLim', [0 255]); %'YLim',[0 0.1*10^4]);


