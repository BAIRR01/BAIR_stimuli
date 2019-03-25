function [newDisplay , selectionMade] = bairUseAlternateDisplay(chosenDisplay)
% Checks if user would like to use an alternate display for experiment
prompt = sprintf('Current display : \n%s \n\nWould you like to use an alternate display? Y/N ', chosenDisplay);
answer = inputdlg(prompt,'Input', 1,{'n'});
if string(answer) == 'Y'|| string(answer) == 'y'
    displays    = {'BAIR_ACER'; 'UMC_OR'};
    
    [whichDisplay, selectionMade] = listdlg('PromptString', 'Which display?', 'SelectionMode', 'single', 'ListString', displays);
    newDisplay = displays{whichDisplay};
else 
    newDisplay = [];
    selectionMade = false;
end
