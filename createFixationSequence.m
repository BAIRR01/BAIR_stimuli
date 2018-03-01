function fixSeq = createFixationSequence(stimulus, dwellTimePerImage, minDurationInSeconds, maxDurationInSeconds)

% Create the fixation dot color change sequence
fixSeq       = ones(size(stimulus.seqtiming));

this_frame = 0;

minDurationInImageNumber = round(minDurationInSeconds / dwellTimePerImage);
maxDurationInImageNumber = round(maxDurationInSeconds / dwellTimePerImage);

while true
    % wait between minDurationInSeconds and maxDurationInSeconds before
    % flipping the dot color
    isi = randperm(maxDurationInImageNumber-minDurationInImageNumber,1)+minDurationInImageNumber-1;
    this_frame = this_frame + isi;
    if this_frame > length(fixSeq), break; end
        fixSeq(this_frame:end) = mod(fixSeq(this_frame-1),2)+1;
end
