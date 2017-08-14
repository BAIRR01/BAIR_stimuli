function BAIR_MEG(n, stimfile)
%RUNME BAIR_MEG(n, stimfile)
%
% MEG BAIR spatiotemporal experiments
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
%   BAIR_MEG(1, 'spatiotemporal_MEG_');


%% 

if notDefined('n'), n = 1; end
if notDefined('stimfile'), stimfile = 'spatiotemporal_MEG_'; end

% initialize stim tracker for MEG
PTBInitStimTracker;
global PTBTriggerLength 
PTBTriggerLength = 0.001;

% debug mode?
% PsychDebugWindowConfiguration
Screen('Preference', 'SkipSyncTests', 1);

%% Calibration
cal = 'meg_lcd';
d   = loadDisplayParams(cal);

% Do we want to use the eyetracker?
use_eyetracker = false;

% d = openScreen(d);
% global PTBTheWindowPtr
% PTBTheWindowPtr = d.windowPtr;

if use_eyetracker
    
    %Open the screen

    
        PTBInitEyeTracker;
%         paragraph = {'Eyetracker initialized.','Get ready to calibrate.'};
%         PTBDisplayParagraph(paragraph, {'center',30}, {'a'});
        PTBCalibrateEyeTracker;
        
        % actually starts the recording
        % name correponding to MEG file (can only be 8 characters!!, no extension)
        PTBStartEyeTrackerRecording('eyelink');
end

%% Default parameters
params = retCreateDefaultGUIParams;


%% Hemifield and ONOFF mixture
params.modality         = 'MEG'; 
params.prescanDuration  = 0;
params.calibration      = cal;
params.startScan        = 0;
params.repetitions      = 1;
params.experiment       = 'Experiment From File';
params.saveMatrix       = 'saveMe';
params.skipSyncTests    = false;

switch stimfile
    case {'spatiotemporal_MEG_' 'ret_MEG_' 'hrf_MEG_'}
        params.fixation = 'disk';
    case 'task_MEG_'
        params.fixation = '4 color dot';
    otherwise 
        error('Unknown stimfile: %s.', stimfile);
end


%% ********************
%  ***** GO ***********
%  *********************
params.loadMatrix = sprintf('%s%d.mat', stimfile, n);
ret(params);

if use_eyetracker

    % retrieve the file
    %     PTBDisplayParagraph({'The experiment is now over.','Please lie still while we save the data.'}, {'center', 30}, {.1});
    
    % Not for now
    %     PTBStopEyeTrackerRecording; % <----------- (can take a while)
    %
    %     % move the file to the logs directory
    %     destination = [pwd '_eyelink_'];
    %     i = 0;
    %     while exist([destination num2str(i) '.edf'], 'file')
    %         i = i + 1;
    %     end
    %     movefile('eyelink.edf', [destination num2str(i) '.edf'])

end

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