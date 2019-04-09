
% Which site?
[experimentSpecs] = bairExperimentSpecs('prompt', false);
whichSite = 7;

% Which subject and session?
subjID = 'intraop_P009';
sessionID = '01';

% Which experiments to run?
experimentTypes = {'prf', 'temporalpattern', 'temporalpattern',  'temporalpattern', 'temporalpattern',  'temporalpattern', 'temporalpattern'};
numberOfExperiments = length(experimentTypes);
runIDs = [1 1 2 1 2 1 2];

% Site-specific stuff to do before starting experiment?
checkforSiteSpecificRequest(experimentSpecs,whichSite);

% Run these experiments!
for ii = 1:numberOfExperiments
    quitProg = BAIR_RUNME(lower(experimentTypes{ii}), runIDs(ii), experimentSpecs(whichSite,:), subjID, sessionID);
    if quitProg, break; end
end