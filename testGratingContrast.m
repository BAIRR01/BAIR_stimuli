PAT = createPatternStimulus(stimParams, 30, 1);
PAT8b = uint8((PAT+.5)*255);

im_cell = cell([3 1]);

cyclesPerDegree     = 3;    
peak2peakContrast   = 60;   

% luminance values need to range between -0.5 and 0.5 prior to 8Bit conversion 
for cc = 1:3
    if cc == 1
        gratingOrientation = pi/2; % GRATING: horizontal
    elseif cc == 2
        gratingOrientation = [pi/2 pi]; % PLAID: horizontal + vertical
    elseif cc == 3
        gratingOrientation = (1:16)/pi; % CIRCULAR: 16 superimposed
    end
    
    imageForThisTrial = createGratingStimulus(stimParams,cyclesPerDegree,gratingOrientation, peak2peakContrast);
    im_cell{cc} = imageForThisTrial;
end

GRAT = im_cell{1};
GRAT8b = uint8((GRAT+.5)*255);
PLAID = im_cell{2};
PLAID8b = uint8((PLAID+.5)*255);
CIRC = im_cell{3};
CIRC8b = uint8((CIRC+.5)*255);

figure;hold on

subplot(4,3,1);
imshow(PAT8b); title('highest density pattern');
subplot(4,3,2); 
histogram(PAT);
set(gca, 'XLim', [-1.2 1.2]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(PAT(:)), peak2peak(PAT(:))));
subplot(4,3,3); 
histogram(PAT8b);
set(gca, 'XLim', [0 255]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(PAT8b(:)), peak2peak(PAT8b(:))));

subplot(4,3,4);
imshow(GRAT8b); 
title(sprintf('Michelson %d percent', peak2peakContrast))
subplot(4,3,5); 
histogram(GRAT);
set(gca, 'XLim', [-1.2 1.2]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(GRAT(:)), peak2peak(GRAT(:))));
subplot(4,3,6); 
histogram(GRAT8b);
set(gca, 'XLim', [0 255]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(GRAT8b(:)), peak2peak(GRAT8b(:))));

subplot(4,3,7);
imshow(PLAID8b); 
title(sprintf('Michelson %d percent', peak2peakContrast))
subplot(4,3,8); 
histogram(PLAID);
set(gca, 'XLim', [-1.2 1.2]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(PLAID(:)), peak2peak(PLAID(:))));
subplot(4,3,9); 
histogram(PLAID8b);
set(gca, 'XLim', [0 255]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(PLAID8b(:)), peak2peak(PLAID8b(:))));

subplot(4,3,10);
imshow(CIRC8b); 
title(sprintf('Michelson %d percent', peak2peakContrast))
subplot(4,3,11); 
histogram(CIRC);
set(gca, 'XLim', [-1.2 1.2]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(CIRC(:)), peak2peak(CIRC(:))));
subplot(4,3,12); 
histogram(CIRC8b);
set(gca, 'XLim', [0 255]);
title(sprintf('rms = %.2f peak2peak = %.2f', rms(CIRC8b(:)), peak2peak(CIRC8b(:))));

