function BAIR_RUNME(runNumber, stimPrefix, siteSpecs, subjID)
% BAIR_RUNME(n, stimPrefix, specs, subjID)
%
% Run BAIR experiments (Do not call this function directly. It gets called
% from the wrapper function, s_runme_BAIR)
% ------
%   Run time per experiment = XX seconds
%
% INPUTS
%   runNumber       run number
%   stimPrefix      prefix for the stimulus files containing images
%                      should be 
%                       - spatiotemporal
%                       - task
%                       - hrfchecker
%                       etc
%   siteSpecs           one-row table generated from the function bairExperimentSpecs
%   subjID          alphanumeric subject ID 
%
%   Example
%    experimentSpecs = bairExperimentSpecs;
%    siteSpecs = experimentSpecs(2,:);
%    runnum = 1;
%    stimPrefix = 'hrfchecker';
%    subjID     = 'wl001';
%    BAIR_RUNME(runnum, stimPrefix, siteSpecs, subjID);

if notDefined('runNumber'), runNumber = 1; end
if notDefined('stimPrefix')
    help(mfilename)
    error('stimPrefix is a required input');
end
if notDefined('siteSpecs')
    help(mfilename)
    error('siteSpecs is a required input');
end

% params = retCreateDefaultGUIParams();

% Set parameters for this experiment
params.experiment       = stimPrefix;
params.subjID           = subjID;
params.loadMatrix       = sprintf('%s_%s_%d.mat', stimPrefix, siteSpecs.sites{1}, runNumber);
params.modality         = siteSpecs.modalities{1};
params.site             = siteSpecs.sites{1};
params.calibration      = siteSpecs.displays{1};
params.triggerKey       = siteSpecs.trigger{1};
params.useSerialPort    = siteSpecs.serialport{1};

% Additional parameters 
params.prescanDuration  = 0;
params.startScan        = 0;
params.stimSize        = 'max';
%params.countdown       = 0;
%params.startScan       = 0;
%params.trigger         = 'Scanner triggers computer';
%params.savestimparams  = 1;
%params.repetitions     = 1;
params.runPriority     = 7;

% Specify task for subject
if contains(stimPrefix, 'task') 
    params.fixation = '4 color dot';
else
    params.fixation = 'disk';
end

% Debug mode?
params.skipSyncTests = 1;

% Go!
doExperiment(params);

end
