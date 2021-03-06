function [experimentSpecs, whichSite, selectionMade] = bairExperimentSpecs(varargin)
% Choose a display for running experiments, making stimuli, analyzing data
% 
% Optional input is paired value, {'prompt', [true or false]}
% If prompt is true, then ask the user which display to use

% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('prompt',false, @islogical);
p.parse(varargin{:});
prompt = p.Results.prompt;

% These are the available displays
sites       = {'nyu3t'; 'nyumeg'; 'nyuecog'; 'umc3t'; 'umc7t'; 'umcecog'; 'umcor'; 'master'};
displays    = {'CBI_Propixx'; 'meg_lcd'; 'SoMMacBook'; 'UMC_3TLCD'; 'UMC_7TDLP'; 'UMC_EMU'; 'UMC_OR'; 'HiResDefault'};
modalities  = {'fMRI'; 'MEG'; 'ECoG'; 'fMRI'; 'fMRI'; 'ECoG'; 'ECoG'; 'none'};
radii       = [12.4; 11; 11.8; 8.3; 6.4287; 11.8; 9.2; 8.3];
trigger     = {'5'; '5'; '5'; 49; 49; '5'; '5'; ''};
serialport  = [false; false; false; true; true; true; true; false];
eyetracker  = [false; true; false; false; false; false; false; false];
photodiode  = [false; false; true; false; false; false; false; false];
displaypos  = [[0 0]; [0 0];[0 0];[0 0];[0 100];[0 0];[0 0]; [0 0]];

experimentSpecs = table(displays, ...
    modalities, ... 
    radii, ...          degrees of visual angle from center of screen to top
    trigger, ...        trigger key to start experiment
    serialport, ...     boolean
    eyetracker, ...     boolean
    photodiode, ...     boolean
    displaypos, ...     [horizontalShift verticalShift] from screen center in pixels
    sites, ...
    'RowNames', sites);

% If requested, ask the user to select a display
if prompt
    [whichSite, selectionMade] = listdlg('PromptString', 'Which site?', 'SelectionMode', 'single', 'ListString', sites);
end

end
