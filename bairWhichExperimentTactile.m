function [experimentType, selectionMade] = bairWhichExperimentTactile

% Which experiment?
experimentTypes = {'SIMPLEHANDSWEEP', 'BLOCKED', 'TACTILEVISUALSWEEP', 'TACTILEVISUALBLOCKED' };
% 'HRFCHECKER'  'HRFPATTERNINVERTED' 'HRFCHECKERINVERTED' };
[whichExperiment, selectionMade] = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);

if selectionMade 
    experimentType = experimentTypes{whichExperiment}; 
else
    experimentType = '';
end
   
end