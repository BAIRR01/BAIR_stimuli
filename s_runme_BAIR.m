
% Which site?
[experimentSpecs, whichSite, ok] = bairExperimentSpecs('prompt', true);
if ~ok, return; end

siteSpecs = experimentSpecs(whichSite,:);

% Which experiment to run?
[experimentType, numberOfRuns] = bairWhichExperiment();

% Prompt for patient ID
prompt = {'Enter subject ID'};
defaults = {'Test099'};
[answer, ok] = inputdlg(prompt, 'Subject Number', 1, defaults);
if ~ok, return; end

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
for n = 1:numberOfRuns
    BAIR_RUNME(n, lower(experimentType), siteSpecs, subjID)    
end






