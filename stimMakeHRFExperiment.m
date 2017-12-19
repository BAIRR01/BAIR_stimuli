function stimMakeHRFExperiment(stimParams, runnum, stimDurationSeconds, repeatStyle)
%stimMakeHRFExperiment(stimParams)
%  Stimuli are presented for 0.25 seconds, with an exponentially
%  distributed interstimulus interval (mean ~9s, range 3-24s).

%% Experiment timing
% stimDurationSeconds = 0.125; % seconds
% repeatStyles = {'same' 'inverted' 'random' 'none'};



numStimuli  = 32;
minISI      = 3;
maxISI      = 24;

% Choose ISIs with an exponential distribution, ranging from minISI to
% maxISI in numStimuli equal steps
[onsets, indices] = getExponentialOnsets(numStimuli, minISI, maxISI, stimDurationSeconds);

%% Spatial patterns

relCutoff = 1/20;
sz = size(stimParams.stimulus.images);
numImages = numStimuli*2+1;

BLANK     = numImages;

bkgrnd = 128;
images = zeros([sz numImages], 'uint8')+bkgrnd;
for ii = 1:numStimuli
    
    if contains(repeatStyle, 'checker')
        im = createCheckerboard(stimParams, 1/relCutoff/2);
    else
        im = createPatternStimulus(stimParams, relCutoff);
    end
    
    mask = mkDisc(size(im), size(im,1)/2);
    
    im = im .*mask;
    images(:,:,ii) = uint8((im+.5)*255);
    images(:,:,ii+numStimuli) = uint8((-im+.5)*255);
end



stimulus = [];
stimulus.cmap         = stimParams.stimulus.cmap;
stimulus.srcRect      = stimParams.stimulus.srcRect;
stimulus.dstRect      = stimParams.stimulus.destRect;
stimulus.images       = images;
stimulus.seqtiming    = 0:stimDurationSeconds:300;
stimulus.seq          = zeros(size(stimulus.seqtiming))+BLANK;

stimOrder = randperm(numStimuli);
stimulus.seq(indices) = stimOrder;

switch repeatStyle
    % repeatStyle = {'same' 'inverted' 'random' 'none'};
    case {'same' 'checkersame'}
        stimulus.seq(indices+1) = stimOrder;  
    case {'none' 'checkernone'}
    % do nothing
    case {'inverted' 'checkerinverted'}
        stimulus.seq(indices+1) = stimOrder+numStimuli;
    case {'random' 'checkerrandom'}
        stimulus.seq(indices+1) = randperm(numStimuli);
    otherwise
        error('Unknown stimulus repeat style %s', repeatStyle);
    
end

stimulus.fixSeq       = ones(size(stimulus.seqtiming));

this_frame = 0;
while true
    % wait between 4 and 20 frames (1 to 5 seconds)
    isi = randperm(16,1)+3;
    this_frame = this_frame + isi;
    if this_frame > length(stimulus.fixSeq), break; end
    stimulus.fixSeq(this_frame:end) = mod(stimulus.fixSeq(this_frame-1),2)+1;
end

switch lower(stimParams.modality{1})
    case 'fmri'
    otherwise
        stimulus.trigSeq  = double(stimulus.seq>0);
        stimulus.diodeSeq = stimulus.trigSeq;
end



fname = sprintf('hrf%s_%s_%d', repeatStyle, stimParams.modality{1}, runnum);
save(fullfile(vistadispRootPath, 'Retinotopy', 'storedImagesMatrices',  fname), 'stimulus')
%plot(stimulus.seqtiming, stimulus.seq)


end

function [onsets, indices] = getExponentialOnsets(numStimuli, minISI, maxISI,temporalResolution)

% Draw numStimuli ISIs from an exponential distribution from [minISI maxISI]
x = linspace(0,1,numStimuli-1) * (1-exp(-minISI)) + exp(-minISI);
ISIs = -log(x)/minISI;
ISIs = ISIs*(maxISI-minISI)+minISI;

% % DEBUG
% disp(mean(ISIs))
%
% figure(1),clf, set(gcf, 'Color', 'w')
% set(gca, 'FontSize', 24); hold on
% plot(ISIs, 'o-', 'LineWidth', 4, 'MarkerSize', 12); axis tight
% ylabel('ITI (s)'); xlabel('Trial')

% Round off the ISIs to multiples of temporalResolution
ISIs = round(ISIs/temporalResolution)*temporalResolution;

% Compute the cummulative sum of ISIs to get the onset times
onsets = cumsum([8 ISIs(randperm(numStimuli-1))]);
indices  = round(onsets/temporalResolution);
% % Debug
% figure(2), clf; set(gcf, 'Color', 'w')
% set(gca, 'FontSize', 24, 'XTick', 0:60:300, 'YTick', []); hold on
% stem(onsets, ones(1,numStimuli+1), 'LineWidth', 2)
% xlabel('Time (s)')
% hgexport(gcf, fullfile(BAIRRootPath, 'figures', 'HRF_onsets.eps'))

end
