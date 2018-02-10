
% Which site?
[experimentSpecs, whichSite] = bairExperimentSpecs('prompt', true);
siteSpecs = experimentSpecs(whichSite,:);

% Prompt for patient ID
prompt = {'Enter subject ID'};
defaults = {'Test099'};
answer = inputdlg(prompt, 'Subject Number', 1, defaults);
subjID = answer{1,:};

% Site-specific stuff
switch siteSpecs.Row{1}
    case 'NYU-ECOG'
        % calibrate display
        NYU_ECOG_Cal();
        
        % Check paths
        if isempty(which('PsychtoolboxVersion'))
            error('Please add Psychtoolbox to path before running')
        end
    otherwise
        % for now, do nothing
end

% Do it!

% spatiotemporal: 176 s for fMRI, XX for ECOG
for n = 101:108
    BAIR_RUNME(n, 'spatiotemporal', siteSpecs, subjID)
end

return

% hrf: 300 SECONDS
for n = 1:1
    BAIR_RUNME(n, 'hrfpattern', siteSpecs, subjID)
    %BAIR_RUNME(n, 'hrfpatterninverted', siteSpecs, subjID)
    %BAIR_RUNME(n, 'hrfcheckerinverted', siteSpecs, subjID)
    %BAIR_RUNME(n, 'hrfchecker', siteSpecs, subjID)
end


% task: 216 SECONDS
for n = 1:10
    BAIR_RUNME(mod(n,2)+1, 'task', siteSpecs, subjID)
end


% retinotopy: 204 SECONDS 
for n = 1:2
    BAIR_RUNME(n, 'ret', siteSpecs, subjID)
end






