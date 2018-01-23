function [experimentSpecs, whichSite] = bairExperimentSpecs(varargin)
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
sites       = {'Master'; 'NYU-3T'; 'NYU-MEG'; 'NYU-ECoG'; 'UMC-3T'; 'UMC-7T'; 'UMC-ECoG'};
displays    = {'HiResDefault'; 'CBI_Propixx'; 'meg_lcd'; 'SoMMacBook'; 'UMC_3TLCD'; 'UMC_7TDLP'; 'default'};
modalities  = {'fMRI'; 'fMRI'; 'MEG'; 'ECoG'; 'fMRI'; 'fMRI'; 'ECoG'};
radii       = {8.3; 12.4; 11; 11.8; 8.3; 6.4287; 11.8};

experimentSpecs = table(displays, modalities, radii, 'RowNames', sites);

% If requested, ask the user to select a display
if prompt
    whichSite = listdlg('PromptString', 'Which site?', 'SelectionMode', 'single', 'ListString', sites);
end

end
