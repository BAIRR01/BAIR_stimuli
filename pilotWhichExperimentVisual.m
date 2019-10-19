function [experimentType, selectionMade] = pilotWhichExperimentVisual

% Which experiment?
experimentTypes = {'SIXCATLOC','SIXCATLOCTEMPORAL','OBJECTDETECTION'}; 
%'EIGHTCATLOC' 'FLOC' 'SCENEFACELATERAL' 'KRAVITZSCENES','BONNERSCENES'
[whichExperiment, selectionMade] = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);

if selectionMade 
    experimentType = experimentTypes{whichExperiment}; 
else
    experimentType = '';
end
   
end