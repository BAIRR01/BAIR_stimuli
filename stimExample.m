function s_example = stimExample(scrnsize)
% Get the data structure for an example stimulus. The fields will be
% populated for actual stimuli.
%
% s_example = stimExample([1024 768]);
%

s_height = min(scrnsize);
s_width  = max(scrnsize);
imsize   = s_height;

% colormap: 1:256 (gamma table is controlled separately)
stimulus.cmap    = (ones(3,1)*(1:256))'; 

% image size and position
stimulus.srcRect = [0 0 imsize imsize];
stimulus.dstRect = [(s_width-imsize)/2   0  (s_width+imsize)/2  imsize];

% a blank image
stimulus.images = uint8(zeros(imsize)+127);

stimulus.seqtiming  = [];
stimulus.seq        = [];
stimulus.fixSeq     = [];


s_example.stimulus = stimulus;

return

% % load example file, stored on the web
% readPth  = 'https://wikis.nyu.edu/download/attachments/85394548/BAIRstimExample.mat?api=v2';
% 
% % local directory where it will be copied
% stimDir  = fullfile(BAIRRootPath, 'stimuli');
% fname    = 'BAIRstimExample.mat';
% writePth = fullfile(stimDir, fname);
% 
% % check whether local directory exists
% if ~exist(stimDir, 'dir'), mkdir(stimDir); end
% addpath(stimDir);
% 
% % copy from web to local
% websave(writePth,readPth);
% 
% % load it into the workspace
% s_example = load(writePth);
