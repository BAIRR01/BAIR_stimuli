% Make Visual Experiment Files
%
% 1. HRF
%       Exponential ISI, mean 9, range 3-24 s
% 2. TASK
%       12 s fixation task, 12 s relax, repeat 8 times (to measure anticipatory BOLD)
% 3. RETINOTOPY
%       Like Dumoulin and Wandell, 2008: 8 sweeps: 4 cardinal, 4 diagonal (diagonals include 50% blanks)
% 4. SPATIOTEMPORAL
%        VISUAL: 36 unique stimuli, shown once each per scan (0.5 s except for temporal stimuli),
%                with mean ISI of 4.5 s, range 3-6 s; orientation (3; 1 grating, 1 plaid, 1 circular);
%                contrast (5; noise patterns); spacing: (5: noise patterns, 1 overlaps with contrast);
%                objects (12: 3 faces, 3 scenes, 3 objects, 3 bodies); temporal (12; 6 durations; 6 ISIs);


% QUESTIONS FOR UTRECHT
%
% HOW DO YOU DEAL WITH TRIGGERS (FMRI)/ ECOG?
% CALIBRATION FILES (INC SCREEN RESOLUTION)
% HOW DO YOU RECEIVE KEY PRESSES? (65-68 for fMRI)
% CHECK ON GAMMA TABLES INC NYU SOM, UTRECHT IEMU, 3T/7T

%   max stimulus radius (in deg)
%       16.6º is the height of the screen for the 3T display at Utrecht,
%       which is the smallest FOV among NYU-3T, UMC-3T, NYU-ECoG, UMC-ECoG,
%       NYU-MEG

stimDiameterDeg = 16.6;       % degrees
peakSFcpd       = 3;          % peak sf of all stimuli (and therefore peak of bandpass filter to make stimuli)
sfAtHalfMax     = [1.4 4.7];  % spatial frequencies where filter falls off to half-height

numberOfRuns    = 1;         

[experimentSpecs, whichSite] = bairExperimentSpecs('prompt', true);

% Generate stimulus template
stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);

stimParams.bpFilter = stimMakeBandPassFilter(stimParams, peakSFcpd, sfAtHalfMax);

% Directory for storing stim files
stimdir = fullfile(BAIRRootPath, 'stimuli');
if ~exist(stimdir, 'dir'), mkdir(stimdir); end

% Make HRF experiment
%   Timing should be specified with values that are integer multiples of
%   the default refresh rate of 60 Hz (ie 16.66666 ms). And the stimulus
%   duration should be an even multiple of this value since we may use
%   paired stimuli for the hRF experiment (a contrast pattern immediately
%   followed by its contrast-reversed pattern). 
stimulusDuration  = 0.200; % seconds. 
dwellTimePerImage = 0.050; % temporal resolution, in s, at which the image sequence is specified 

for runNum = 1:numberOfRuns
    stimMakeHRFExperiment(stimParams, runNum, stimulusDuration, dwellTimePerImage,  'pattern');
%     stimMakeHRFExperiment(stimParams, runNum, stimulusDuration, dwellTimePerImage,  'patternInverted');
%     stimMakeHRFExperiment(stimParams, runNum, stimulusDuration, dwellTimePerImage,  'checker');
%     stimMakeHRFExperiment(stimParams, runNum, stimulusDuration, dwellTimePerImage,  'checkerinverted');
end

return

% MAKE TASK EXPERIMENT
stimMakeTaskExperiment(stimParams, 'fMRI');
stimMakeTaskExperiment(stimParams, 'MEG');

% Make PRF experiment
stimMakePRFExperiment(stimParams, 'fMRI');
stimMakePRFExperiment(stimParams, 'MEG');

% Make spatiotemporal experiment
stimMakeSpatiotemporalExperiment(stimParams, 'fMRI');
stimMakeSpatiotemporalExperiment(stimParams, 'MEG');


