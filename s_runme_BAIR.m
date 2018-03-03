
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

% % Which experiment to run?
% [experimentType] = bairWhichExperiment();
% [runID] = bairWhichRun();

% Which experiments to run?
fname = 'experimentsToRun';
T = readtable(fullfile(vistadispRootPath, sprintf('%s.txt', fname)));

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

% Run these experiments!
numberOfExperiments = height(T);
for ii = 1:numberOfExperiments
    experimentType = T.Var1{ii};
    runID = T.Var2(ii);
    BAIR_RUNME(lower(experimentType), runID, siteSpecs, subjID)    
end






