function [experimentType, selectionMade] = bairWhichExperiment

% Which experiment?
experimentTypes = {'DOTTASK' 'PRF' 'SPATIALPATTERN' 'SPATIALOBJECT' ...
    'TEMPORALPATTERN' 'HRFPATTERN' 'HRFPATTERNBREATHINGCHALLENGE','TACTILE',...
    'FINGERMAPPINGLEFT','FINGERMAPPINGRIGHT' 'GESTURES', 'GESTURESLEARNING', ...
    'GESTURESPRACTICE','BOLDHAND', 'BOLDSAT' };
% 'HRFCHECKER'  'HRFPATTERNINVERTED' 'HRFCHECKERINVERTED' };
[whichExperiment, selectionMade] = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);

if selectionMade 
    experimentType = experimentTypes{whichExperiment}; 
else
    experimentType = '';
end
   
end
