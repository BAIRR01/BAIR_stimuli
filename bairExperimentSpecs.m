function [experimentSpecs, whichSite, ok] = bairExperimentSpecs(varargin)
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
sites       = {'NYU-3T'; 'NYU-MEG'; 'NYU-ECOG'; 'UMC-3T'; 'UMC-7T'; 'UMC-ECOG';'Master'};
displays    = {'CBI_Propixx'; 'meg_lcd'; 'SoMMacBook'; 'UMC_3TLCD'; 'UMC_7TDLP'; 'default'; 'HiResDefault'};
modalities  = {'fMRI'; 'MEG'; 'ECoG'; 'fMRI'; 'fMRI'; 'ECoG'; 'none'};
radii       = {12.4; 11; 11.8; 8.3; 6.4287; 11.8; 8.3};
trigger     = {'5'; '5'; '5'; 49; 49; '5'; ''};
serialport  = {false; false; false; true; true; true; false};
eyetracker  = {false; true; false; false; false; false; false};
    
experimentSpecs = table(displays, modalities, radii, trigger, serialport, eyetracker, sites, 'RowNames', sites);

% If requested, ask the user to select a display
if prompt
    [whichSite, ok] = listdlg('PromptString', 'Which site?', 'SelectionMode', 'single', 'ListString', sites);
end

end
