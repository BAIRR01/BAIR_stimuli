function [experimentType, selectionMade] = pilotWhichExperimentVisual

% Which experiment?
experimentTypes = {'SIXCATLOC' 'EIGHTCATLOC' 'KALANITLOC' ...
    'SCENEFACELATERAL' 'KRAVITZSCENES','BONNERSCENES', ...
    'OBJECTDETECTION'};
[whichExperiment, selectionMade] = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);

if selectionMade 
    experimentType = experimentTypes{whichExperiment}; 
else
    experimentType = '';
end
   
end