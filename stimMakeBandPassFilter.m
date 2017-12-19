function bandpassFilter   = stimMakeBandPassFilter(stimParams, peakSFcpd, sfAtHalfMax)
% Make bandpass filter for generating BAIR visual stimuli
%
% Filter properties are based on Kay et al, 2013:
% "Stimuli consisted of grayscale images restricted to a band-pass
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
% naturalistic stimuli."


stimDiameterInPixels = stimParams.stimulus.srcRect(4);

% % Make a bandpass filter from a laplacian of gaussian
% %bpfilter = -fspecial('log', 50, 4);
% 
% Make a bandpass filter from a difference of two gaussians
%       The constants 0.011 and 1.5 were determined by trial and error to
%       give the approximately correct results
centerSigma     = 0.011 / peakSFcpd * stimDiameterInPixels;
surroundSigma   = centerSigma * 1.5;
gaussianSupport = round(10 * surroundSigma);


centerFilter   = fspecial('gaussian', gaussianSupport, centerSigma); % high frequency cut-off
surroundFilter = fspecial('gaussian', gaussianSupport, surroundSigma);  % low  frequency cut-off

bandpassFilter =  centerFilter/(mean(centerFilter(:))) - surroundFilter/mean(surroundFilter(:));
bandpassFilter = bandpassFilter / norm(bandpassFilter(:));

% DEBUG
% % Check results with a plot
% figure
% subplot(1,2,1)
% surf(bandpassFilter)
% 
% subplot(1,2,2); hold on
% 
% [frequencies, amplitudes, binnedFrequencies, binnedAmplitudes] = ...
%     compute2DamplitudeSpectrum(bandpassFilter, stimParams.display);
% 
% plot(frequencies(:), amplitudes(:), '-', binnedFrequencies, binnedAmplitudes);
% set(gca, 'XTick', sort([peakSFcpd sfAtHalfMax]), 'XGrid', 'on', 'XLim', [0 10])
% xlabel('Frequency (cycles per degree)')



