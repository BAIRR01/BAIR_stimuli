% task
% Make sure vistadisp is set to master branch

% 204 SECONDS, 1.5 S TR (136 VOLUMES)
for n = 1:2
    BAIR_FMRI(n, 'ret_fMRI_')
end

% 176 SECONDS, 1.5 S TR (118 VOLUMES)
for n = 1:8
    BAIR_FMRI(n, 'spatiotemporal_fMRI_')
end

% 216 SECONDS, 1.5 S TR (144 VOLUMES)
for n = 1:2
    BAIR_FMRI(n, 'task_fMRI_')
end

% load /Users/winawerlab/matlab/git/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices/ret_fMRI_2.mat
% idx = sort(randperm(408, 50));
% for ii = 1:50
%     stimulus.fixSeq(idx(ii):end) = mod(ii,2)+1;
% end
% 
% save /Users/winawerlab/matlab/git/vistadisp/Applications2/Retinotopy/standard/storedImagesMatrices/ret_fMRI_2.mat stimulus




