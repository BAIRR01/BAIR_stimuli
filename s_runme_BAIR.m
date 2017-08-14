
modalities = {'fMRI' 'MEG' 'EEG' 'ECoG'};
m = input('Which modality (1/2/3/4)? 1 fMRI; 2 MEG; 3 EEG; 4 ECoG\n');
modality = modalities{m};

switch modality
    case 'fMRI', runmefun = @(n, str) BAIR_FMRI(n, str);
    case 'MEG',  runmefun = @(n, str) BAIR_MEG(n, str);
    otherwise, error('Not yet implemented');
end

runmefun(n, sprintf('hrf_%s_', modality))
runmefun(mod(n,2)+1, sprintf('task_%s_', modality))
runmefun(n, sprintf('ret_%s_', modality))
runmefun(n, sprintf('spatiotemporal_%s_', modality))

% hrf: 300 SECONDS
for n = 1:12
    runmefun(n, sprintf('hrf_%s_', modality))
end

% 216 SECONDS
for n = 1:10
    runmefun(mod(n,2)+1, sprintf('task_%s_', modality))
end


% retinotopy: 204 SECONDS 
for n = 1:2
    runmefun(n, sprintf('ret_%s_', modality))
end

% 176 SECONDS
for n = 1:8
    runmefun(n, sprintf('spatiotemporal_%s_', modality))
end





