
% Which site?
[experimentSpecs, whichSite] = bairExperimentSpecs('prompt', true);
siteSpecs = experimentSpecs(whichSite,:);

% Prompt for patient ID
prompt = {'Enter subject ID'};
defaults = {'Test099'};
answer = inputdlg(prompt, 'Subject Number', 1, defaults);
subjID = answer{1,:};


% Name of modality-specific runme function
switch lower(siteSpecs.modalities{1})
    % a: run number (integer)
    % b: stimulus prefix
    % c: site specs
    % d: subject ID
    case 'fmri', runmefun = @(a,b,c,d) BAIR_FMRI(a,b,c,d);
    case 'meg',  runmefun = @(a,b,c,d) BAIR_MEG(a,b,c,d);
    case 'ecog', runmefun = @(a,b,c,d) BAIR_ECOG(a,b,c,d);        
    otherwise, error('Not yet implemented');
end

% Do it!

% hrf: 300 SECONDS
for n = 1:1
    runmefun(n, 'hrfpatterninverted', siteSpecs, subjID)
    runmefun(n, 'hrfpattern', siteSpecs, subjID)
    runmefun(n, 'hrfcheckerinverted', siteSpecs, subjID)
    runmefun(n, 'hrfchecker', siteSpecs, subjID)
end

% RUN ONLY HRF for now
return

% spatiotemporal: 176 SECONDS
for n = 1:8
    runmefun(n, 'spatiotemporal', siteSpecs, subjID)
end

% task: 216 SECONDS
for n = 1:10
    runmefun(mod(n,2)+1, 'task', siteSpecs, subjID)
end


% retinotopy: 204 SECONDS 
for n = 1:2
    runmefun(n, 'ret', siteSpecs, subjID)
end






