function stimulus = sparsifyStimulusStruct(stimulus, maxUpdateInterval)

numImagesIn = length(stimulus.seq);
keepFrames = true(1,numImagesIn);

mostRecentFrame = 1;

if ~isfield(stimulus, 'trigSeq')
    stimulus.trigSeq = zeros(size(stimulus.seq));
    removeTrigSeq = true;
else
    removeTrigSeq = false;
end

for ii = 2:numImagesIn - 1
   thisInterval = stimulus.seqtiming(ii) - stimulus.seqtiming(mostRecentFrame);
      
   didAnyThingChange = any([...
       stimulus.seq(ii) ~= stimulus.seq(mostRecentFrame) ...
       stimulus.fixSeq(ii) ~= stimulus.fixSeq(mostRecentFrame) ...
       stimulus.trigSeq(ii) ~= stimulus.trigSeq(mostRecentFrame) ...
       ]);
   
   if thisInterval < maxUpdateInterval && ~didAnyThingChange
      keepFrames(ii) = false; 
      
   else
       mostRecentFrame = ii;
   end
end
   
stimulus.seqtiming = stimulus.seqtiming(keepFrames);
stimulus.seq       = stimulus.seq(keepFrames);
stimulus.fixSeq    = stimulus.fixSeq(keepFrames);
stimulus.trigSeq   = stimulus.trigSeq(keepFrames);

if removeTrigSeq, stimulus = rmfield(stimulus, 'trigSeq'); end
end