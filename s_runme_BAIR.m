
[experimentSpecs, whichSite] = bairExperimentSpecs('prompt', true);

modality = experimentSpecs.modalities{whichSite};

switch lower(modality)
    case 'fmri', runmefun = @(n, str) BAIR_FMRI(n, str);
    case 'meg',  runmefun = @(n, str) BAIR_MEG(n, str);
    case 'ecog', runmefun = @(n, str) BAIR_ECOG(n, str);        
    otherwise, error('Not yet implemented');
end

for n = 1:2
    runmefun(n, sprintf('hrfpatterninverted_%s_', modality))
    
    return
    runmefun(n, sprintf('hrfpattern_%s_', modality))
    runmefun(n, sprintf('hrfcheckerinverted_%s_', modality))
    runmefun(n, sprintf('hrfchecker_%s_', modality))
end

% hrf: 300 SECONDS

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





