% Which site?
[experimentSpecs, whichSite, selectionMade] = bairExperimentSpecs('prompt', true);
if ~selectionMade, return; end

% Which subject and session?
[subjID, sessionID, ssDefined] = bairWhichSubjectandSession();
if ~ssDefined, return; end

% Which type of experiment?
[sensoryDomain, selectionMade] = bairWhichSensoryModality();
if ~selectionMade, return; end

% Which experiments to run?
[numberOfExperiments, experimentTypes, runIDs, fileSelected] = bairWhichExperimentList(experimentSpecs.sites{whichSite});
if ~fileSelected, return; end

% Site-specific stuff to do before starting experiment?
checkforSiteSpecificRequest(experimentSpecs, whichSite, sensoryDomain);

% Run these experiments!
for ii = 1:numberOfExperiments
    quitProg = BAIR_RUNME(lower(experimentTypes{ii}), runIDs(ii), experimentSpecs(whichSite,:), subjID, sessionID, sensoryDomain);
    if quitProg, break; end
end