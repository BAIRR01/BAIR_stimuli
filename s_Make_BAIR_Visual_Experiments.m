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


% These are the available displays
sites       = {'NYU-3T'; 'NYU-MEG'; 'NYU-ECoG'; 'UMC-3T'; 'UMC-7T'; 'UMC-ECoG'};
displays    = {'CBI_Propixx'; 'meg_lcd'; 'SoMMacBook'; 'default'; 'default'; 'default'};
modalities  = {'fMRI'; 'MEG'; 'ECoG'; 'fMRI'; 'fMRI'; 'ECoG'};
%radii       = {12.4 11 11.8 8.3 6.45 11.8}

experimentSpecs = table(displays, modalities, 'RowNames', sites);

% USER SPECIFIED
%   which display?
whichSite = listdlg('PromptString', 'Which site?', 'SelectionMode', 'single', 'ListString', sites);

% Generate stimulus template
stimParams = stimInitialize(experimentSpecs, whichSite, stimDiameterDeg);

stimParams.bpFilter = stimMakeBandPassFilter(stimParams, peakSFcpd, sfAtHalfMax);

% Directory for storing stim files
stimdir = fullfile(BAIRRootPath, 'stimuli');
if ~exist(stimdir, 'dir'), mkdir(stimdir); end

% Make HRF experiment
stimMakeHRFExperiment(stimParams, 1, .125, 'same');
stimMakeHRFExperiment(stimParams, 2, .125, 'same');
stimMakeHRFExperiment(stimParams, 3, .125, 'same');
stimMakeHRFExperiment(stimParams, 4, .125, 'same');
stimMakeHRFExperiment(stimParams, 1, .125, 'inverted');
stimMakeHRFExperiment(stimParams, 2, .125, 'inverted');
stimMakeHRFExperiment(stimParams, 3, .125, 'inverted');
stimMakeHRFExperiment(stimParams, 4, .125, 'inverted');

stimMakeHRFExperiment(stimParams, 1, .125, 'checkersame');
stimMakeHRFExperiment(stimParams, 2, .125, 'checkersame');
stimMakeHRFExperiment(stimParams, 3, .125, 'checkersame');
stimMakeHRFExperiment(stimParams, 4, .125, 'checkersame');
stimMakeHRFExperiment(stimParams, 1, .125, 'checkerinverted');
stimMakeHRFExperiment(stimParams, 2, .125, 'checkerinverted');
stimMakeHRFExperiment(stimParams, 3, .125, 'checkerinverted');
stimMakeHRFExperiment(stimParams, 4, .125, 'checkerinverted');


% MAKE TASK EXPERIMENT
stimMakeTaskExperiment(stimParams, 'fMRI');
stimMakeTaskExperiment(stimParams, 'MEG');


% Make PRF experiment
stimMakePRFExperiment(stimParams, 'fMRI');
stimMakePRFExperiment(stimParams, 'MEG');

% Make spatiotemporal experiment
stimMakeSpatiotemporalExperiment(stimParams, 'fMRI');
stimMakeSpatiotemporalExperiment(stimParams, 'MEG');


