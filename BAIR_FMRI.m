function BAIR_FMRI(n, stimPrefix, specs, subjID)
% BAIR_FMRI(n, stimPrefix, specs, subjID)
%
% FMRI BAIR experiments
% ------
%   Run time per experiment = XX seconds
%
% INPUTS
%   n:              run number
%   stimPrefix      prefix for the stimulus files containing images
%                      should be 
%                       - spatiotemporal
%                       - task
%                       - hrfchecker
%                       etc
%   specs           one-row table generated from the function bairExperimentSpecs
%   subjID          alphanumeric subject ID 
%
%   Example
%    experimentSpecs = bairExperimentSpecs;
%    siteSpecs = experimentSpecs(2,:);
%    runnum = 1;
%    stimPrefix = 'hrfchecker';
%    subjID     = 'wl001';
%    BAIR_FMRI(runnum, stimPrefix, siteSpecs, subjID);


%% 

if notDefined('n'), n = 1; end
if notDefined('stimPrefix')
    help(mfilename)
    error('stimPrefix is a required input');
end
if notDefined('specs')
    help(mfilename)
    error('specs is a required input');
end

% debug mode?
% PsychDebugWindowConfiguration
Screen('Preference', 'SkipSyncTests', 1);

% Calibration
cal = specs.displays{1};

% Site
site = specs.Row{1};

% Default parameters
params = retCreateDefaultGUIParams;


% Set parameters for this experiment
params.modality         = 'fMRI'; 
params.prescanDuration  = 0;
params.calibration      = cal;
params.startScan        = 0;
params.repetitions      = 1;
params.experiment       = 'Experiment From File';
params.triggerKey       = '5';
params.devices          = 'External: 2';
params.saveMatrix       = 'saveMe';
params.loadMatrix       = sprintf('%s_%s_%d.mat', stimPrefix, site, n);
params.skipSyncTests    = false;
params.fixation         = 'disk';
params.subjID           = subjID;

if contains(stimPrefix, 'task'), params.fixation = '4 color dot';end

% Go!
ret(params);

end
