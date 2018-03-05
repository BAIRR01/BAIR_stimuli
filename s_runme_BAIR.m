
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
filespec = fullfile(vistadispRootPath, 'RunMe', '*.txt');
[fname, pathname] = uigetfile(filespec, 'Select a text file with the experiment list');

T = readtable(fullfile(pathname, fname));
experimentType      = T.Var1;
runID               = T.Var2;
numberOfExperiments = height(T);

% Check that the experiment files exist
stimPath = fullfile(vistadispRootPath, 'StimFiles');
for ii = 1:numberOfExperiments
    fname = sprintf('%s_%s_%d.mat', siteSpecs.sites{1}, experimentType{ii}, runID(ii));
    
    if ~exist(fullfile(stimPath, fname), 'file')
        error('Requested experiment file %s not found in expected location:\n%s', fname, stimPath);        
    end

end

fprintf('All experiment files ready to go.\n');

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
    quitProg = BAIR_RUNME(lower(experimentType{ii}), runID(ii), siteSpecs, subjID);
    if quitProg, break; end
end