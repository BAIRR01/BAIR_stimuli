function quitProg = BAIR_RUNME(stimPrefix, runID, siteSpecs, subjID, sessionID, sensoryDomain)
% quitProg = BAIR_RUNME(stimPrefix, runNumber, specs, subjID, sessionID)
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
%   sensoryDomain   should be one of the following sensory modalities
%                         - visual
%                         - tactile
%                         - motor
%
%   Example
%    experimentSpecs = bairExperimentSpecs;
%    siteSpecs       = experimentSpecs(2,:);
%    runnum          = 1;
%    stimPrefix      = 'hrfchecker';
%    subjID          = 'wl001';
%    sensoryDomain   = 'visual'
%    BAIR_RUNME(runnum, stimPrefix, siteSpecs, subjID,[], sensoryDomain);

if notDefined('stimPrefix')
    help(mfilename)
    error('stimPrefix is a required input');
end
if notDefined('sensoryDomain')
    help(mfilename)
    error('sensoryDomain is a required input');
end
if notDefined('siteSpecs')
    help(mfilename)
    error('siteSpecs is a required input');
end
if notDefined('runID'), runID = 1; end
if notDefined('sessionID'), sessionID = '01'; end

% Set parameters for this experiment
params.experiment       = stimPrefix;
params.subjID           = subjID;
params.runID            = runID;
params.sessionID        = sessionID;
params.loadMatrix       = sprintf('%s_%s_%d.mat', siteSpecs.sites{1}, stimPrefix, runID);
params.modality         = siteSpecs.modalities{1};
params.site             = siteSpecs.sites{1};
params.calibration      = siteSpecs.displays{1};
params.triggerKey       = siteSpecs.trigger{1};
params.useSerialPort    = siteSpecs.serialport;
params.useEyeTracker    = siteSpecs.eyetracker;
params.shiftDestRect    = siteSpecs.displaypos;
params.sensoryDomain    = sensoryDomain;
params.useDataGlove     = false; 

% Additional parameters 
params.prescanDuration  = 0;
params.startScan        = 0;

% Set priority (depends on operating system)
if ispc
    params.runPriority  = 2;
elseif ismac
    params.runPriority  = 7;
end

% Specify which type of fixation to use
params.fixation = 'cross'; % default
if contains(stimPrefix, 'dottask') 
    params.fixation = '4 color dot';
end
if contains(lower(sensoryDomain), 'motor')
    params.fixation = 'disk';
end

% Sensory modality-specific stuff to do before starting experiment?
[params] = checkforSensoryDomainSpecificRequest(params);

% Debug mode?
params.skipSyncTests = 1;

% Go!
quitProg = doExperiment(params);

end
