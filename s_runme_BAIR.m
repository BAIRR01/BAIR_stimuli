
modalities = {'fMRI' 'MEG' 'EEG' 'ECoG'};
modality = modalities{1};

switch modality
    case 'fMRI', runmefun = @(n, str) BAIR_FMRI(n, str);
    case 'MEG',  runmefun = @(n, str) BAIR_MEG(n, str);
    otherwise, error('Not yet implemented');
end

% 300 SECONDS
for n = 1:2
    runmefun(n, sprintf('hrf_%s_', modality))
end

% 204 SECONDS
for n = 1:2
    runmefun(n, 'ret_%s_', modality))
end

% 176 SECONDS
for n = 1:8
    runmefun(n, 'spatiotemporal_%s_', modality))
end

% 216 SECONDS
for n = 1:2
    runmefun(n, 'task_%s_', modality))
end



