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


% General

% Size of display in pixels
screensz = [1024 768];

% Example stimulus
s_example = stimExample(screensz);

% Directory for storing stim files
stimdir = fullfile(BAIRRootPath, 'stimuli');
if ~exist(stimdir, 'dir'), mkdir(stimdir); end

% MAKE TASK EXPERIMENT

stimMakeTaskExperiment(s_example, 'fMRI');
stimMakeTaskExperiment(s_example, 'MEG');

% Make HRF experiment
stimMakeHRFExperiment(s_example, 'fMRI');
stimMakeHRFExperiment(s_example, 'MEG');

% Make PRF experiment
stimMakePRFExperiment(s_example, 'fMRI');
stimMakePRFExperiment(s_example, 'MEG');

% Make spatiotemporal experiment
stimMakeSpatiotemporalExperiment(s_example, 'fMRI');
stimMakeSpatiotemporalExperiment(s_example, 'MEG');


