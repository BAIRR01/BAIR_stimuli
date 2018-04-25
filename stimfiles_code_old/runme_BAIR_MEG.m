
% 300 SECONDS
for n = 1:2
    BAIR_MEG(n, 'hrf_MEG_')
end

return

% 204 SECONDS
for n = 1:2
    BAIR_MEG(n, 'ret_MEG_')
end

% 58 SECONDS
for n = 1:8
    BAIR_MEG(n, 'spatiotemporal_MEG_')
end

% 216 SECONDS
for n = 1:2
    BAIR_MEG(n, 'task_MEG_')
end

% load /Users/winawerlab/matlab/git/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices/ret_MEG_2.mat
% idx = sort(randperm(408, 50));
% for ii = 1:50
%     stimulus.fixSeq(idx(ii):end) = mod(ii,2)+1;
% end
% 
% save /Users/winawerlab/matlab/git/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices/ret_MEG_2.mat stimulus




