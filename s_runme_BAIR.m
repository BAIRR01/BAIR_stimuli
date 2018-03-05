
% Which site?
[experimentSpecs, whichSite, ok] = bairExperimentSpecs('prompt', true);
if ~ok, return; end

siteSpecs = experimentSpecs(whichSite,:);

% Prompt for patient ID
prompt = {'Enter subject ID'};
defaults = {'test'};
[answer] = inputdlg(prompt, 'Subject Number', 1, defaults);
if isempty(answer), return; end
subjID = answer{1,:};

% Which experiments to run?
[numberOfExperiments, experimentTypes, runIDs] = bairWhichExperimentList(siteSpecs.sites{1});

% Site-specific stuff
switch siteSpecs.sites{1}
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

% Run these experiments!
for ii = 1:numberOfExperiments
    quitProg = BAIR_RUNME(lower(experimentTypes{ii}), runIDs(ii), siteSpecs, subjID);
    if quitProg, break; end
end