
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

% Which experiment to run?
[experimentType] = bairWhichExperiment();
[runID] = bairWhichRun();

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
%for runNumber = 1:numberOfRuns
    BAIR_RUNME(lower(experimentType), runID, siteSpecs, subjID)    
%end






