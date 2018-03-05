function quitProg = BAIR_RUNME(stimPrefix, runID, siteSpecs, subjID)
% quitProg = BAIR_RUNME(stimPrefix, runNumber, specs, subjID)
%
% Run BAIR experiments (Do not call this function directly. It gets called
% from the wrapper function, s_runme_BAIR)
% ------
%   Run time per experiment = XX seconds
%
% INPUTS
%   stimPrefix      prefix for the stimulus files containing images
%                      should be 
%                       - spatiotemporal
%                       - prf
%                       - spatialpattern
%                       - temporalpattern
%                       etc
%   runID           run number
%   siteSpecs       one-row table generated from the function bairExperimentSpecs
%   subjID          alphanumeric subject ID 
%
%   Example
%    experimentSpecs = bairExperimentSpecs;
%    siteSpecs = experimentSpecs(2,:);
%    runnum = 1;
%    stimPrefix = 'hrfchecker';
%    subjID     = 'wl001';
%    BAIR_RUNME(runnum, stimPrefix, siteSpecs, subjID);

if notDefined('stimPrefix')
    help(mfilename)
    error('stimPrefix is a required input');
end
if notDefined('siteSpecs')
    help(mfilename)
    error('siteSpecs is a required input');
end
if notDefined('runID'), runID = 1; end

% Set parameters for this experiment
params.experiment       = stimPrefix;
params.subjID           = subjID;
params.runID            = runID;
params.loadMatrix       = sprintf('%s_%s_%d.mat', siteSpecs.sites{1}, stimPrefix, runID);
params.modality         = siteSpecs.modalities{1};
params.site             = siteSpecs.sites{1};
params.calibration      = siteSpecs.displays{1};
params.triggerKey       = siteSpecs.trigger{1};
params.useSerialPort    = siteSpecs.serialport;
params.useEyeTracker    = siteSpecs.eyetracker;
params.shiftDestRect    = siteSpecs.displaypos;

% Additional parameters 
params.prescanDuration  = 0;
params.startScan        = 0;
params.stimSize        = 'max';
params.runPriority     = 7;

% Specify task for subject
if contains(stimPrefix, 'task') 
    params.fixation = '4 color dot';
else
    %params.fixation = 'disk';
    params.fixation = 'cross';
end

% Debug mode?
params.skipSyncTests = 1;

% Go!
quitProg = doExperiment(params);

end
