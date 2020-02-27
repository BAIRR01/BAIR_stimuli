function [experimentType, selectionMade] = pilotWhichExperimentVisual

% Which experiment?
experimentTypes = {'SIXCATLOC','SIXCATLOCTEMPORAL','SIXCATLOCTEMPORALDIFF','SIXCATLOCISIDIFF', 'OBJECTDETECTION', 'SCENEFACELATERAL'}; 
%'EIGHTCATLOC' 'FLOC', 'KRAVITZSCENES','BONNERSCENES'
[whichExperiment, selectionMade] = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);

if selectionMade 
    experimentType = experimentTypes{whichExperiment}; 
else
    experimentType = '';
end
   
end