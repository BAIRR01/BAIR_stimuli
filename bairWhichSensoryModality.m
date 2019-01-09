function [sensoryDomain, selectionMade] = bairWhichSensoryModality

% Which experiment?
sensoryModalityTypes = {'VISUAL','TACTILE','MOTOR'};

[bairWhichSensoryModality, selectionMade] = listdlg('PromptString', 'Which sensory modality?', 'ListString', sensoryModalityTypes);

if selectionMade 
    sensoryDomain = sensoryModalityTypes{bairWhichSensoryModality}; 
else
    sensoryDomain = '';
end
   
end