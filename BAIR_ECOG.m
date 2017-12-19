function BAIR_ECOG(n, stimfile)
%RUNME BAIR_ECOG(n, stimfile)
%
% ECoG BAIR spatiotemporal experiments
% ------
%   Run time per experiment = XX seconds
%
% INPUTS
%   n is the runnumber [1 24]
%   stimfile is the prefix for the stimulus fils containing images, and
%           should be
%               - spatiotemporal
%           
% The actual stim files have names like
%   spatiotemporal1.mat
%   spatiotemporal2.mat

% Example
%   BAIR_ECOG(1, 'spatiotemporal_ECOG_');


%% 

if notDefined('n'), n = 1; end
if notDefined('stimfile'), stimfile = 'spatiotemporal_ECOG_'; end

% debug mode?
Screen('Preference', 'SkipSyncTests', 1);

%% Calibration
cal = 'gunjou';

%% Default parameters
params = retCreateDefaultGUIParams;


%% Hemifield and ONOFF mixture
params.modality         = 'ECoG'; 
params.prescanDuration  = 0;
params.calibration      = cal;
params.startScan        = 0;
params.repetitions      = 1;
params.experiment       = 'Experiment From File';
params.saveMatrix       = 'saveMe';
params.skipSyncTests    = false;

switch stimfile
    case {'spatiotemporal_ECOG_' 'ret_ECOG_' 'hrf_ECOG_'}
        params.fixation = 'disk';
    case 'task_ECOG_'
        params.fixation = '4 color dot';
    otherwise 
        error('Unknown stimfile: %s.', stimfile);
end


%% ********************
%  ***** GO ***********
%  *********************
params.loadMatrix = sprintf('%s%d.mat', stimfile, n);
ret(params);
