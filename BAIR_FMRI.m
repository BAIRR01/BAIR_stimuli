function BAIR_FMRI(n, stimfile)
% BAIR_FMRI(n, stimfile)
%
% FMRI BAIR experiments
% ------
%   Run time per experiment = XX seconds
%
% INPUTS
%   n is the runnumber [1 24]
%   stimfile is the prefix for the stimulus fils containing images, and
%           should be
%               - spatiotemporal_fMRI_
%               - task_fMRI_
%           
% The actual stim files have names like
%   spatiotemporal_fMRI_1.mat
%   spatiotemporal_fMRI_2.mat

% Example
%   BAIR_FMRI(1, 'spatiotemporal_fMRI_');


%% 

if notDefined('n'), n = 1; end
if notDefined('stimfile'), stimfile = 'spatiotemporal_fMRI_'; end

% debug mode?
% PsychDebugWindowConfiguration
Screen('Preference', 'SkipSyncTests', 1);

%% Calibration
cal = 'CBI_Propixx';


%% Default parameters
params = retCreateDefaultGUIParams;


%% Hemifield and ONOFF mixture
params.modality         = 'fMRI'; 
params.prescanDuration  = 0;
params.calibration      = cal;
params.startScan        = 0;
params.repetitions      = 1;
params.experiment       = 'Experiment From File';
params.triggerKey       = '5';
params.devices          = 'External: 2';
params.saveMatrix       = 'saveMe';
params.skipSyncTests    = false;

switch stimfile
    case 'task_fMRI_'
        params.fixation = '4 color dot';
    otherwise 
        params.fixation = 'disk';
end



%% ********************
%  ***** GO ***********
%  *********************
params.loadMatrix = sprintf('%s%d.mat', stimfile, n);
ret(params);

%% Check timing results
f = dir('~/Desktop/201*.mat');
load(fullfile('~', 'Desktop', f(end).name));
figure(101); clf

% desired inter-stimulus duration
plot(diff(stimulus.seqtiming));

% measured inter-stimulus duration
hold on; plot(diff(response.flip), 'r-'); 

ylim(median(diff(response.flip)) + [-.001 .001])
% frames between stimuli
frames = round(diff(response.flip) / (1/60)); 

% how many interstimulus frames differed from the median?
disp(sum(frames ~= median(frames)))