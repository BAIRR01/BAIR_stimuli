function [experimentType, selectionMade] = bairWhichExperimentTactile

% Which experiment?
experimentTypes = {'SIMPLEHANDSWEEP', 'BLOCKED', 'TACTILEVISUALSWEEP', 'TACTILEVISUALBLOCKED', 'RANDOM_5FINGERS', 'TEMPORAL', 'STIMULUSTEST'};
% 'HRFCHECKER'  'HRFPATTERNINVERTED' 'HRFCHECKERINVERTED' };
[whichExperiment, selectionMade] = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);

if selectionMade 
    experimentType = experimentTypes{whichExperiment}; 
else
    experimentType = '';
end
   
end