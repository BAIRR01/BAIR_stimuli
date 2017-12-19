%% Load the bpfilter and the base images
%load(fullfile(rootpath, 'data', 'input', 'stimuli.mat'), 'bpfilter');

addpath(genpath('~/matlab/git/soccode'))
addpath(genpath('~/matlab/toolboxes/knkutils'))

%load('bpfilter')
load('/Users/jonathanwinawer/matlab/junk/bpfilter.mat')
%% 
[output, edge, thresh, res] = createPatternStimulus([768, 768], 1/40, bpfilter);
%[output, edge, thresh, res] = createPatternStimulus([768, 768], 1/10, bpfilter);

output = output / max(abs(output(:))) + 0.5;
output(output<0) = 0;output(output>1) = 1;

figure(1), colormap gray
subplot(2,2,1)
imagesc(res); axis image;

subplot(2,2,2)
imagesc(thresh); axis image;

subplot(2,2,3)
imagesc(edge);axis image;

subplot(2,2,4)
imagesc(output);axis image;
%
%% BARS
res         = 768;             % native resolution that we construct at
totalfov    = 12;              % total number of degrees for image
cpd         = 3;               % target cycles per degree
radius      = totalfov/2;      % radius of image in degrees
maskFrac    = 11.5/totalfov;   % what fraction of radius for the circular mask to start at

nframes     = 9;               % how many images from each class

cpim        = totalfov*cpd;    % cycles per image that we are aiming for
spacing     = res/cpim;        % pixels to move from one cycle to the next

span        = [-0.5, 0.5];      % dynamic range, to put into 'imshow' etc.

% Choose which bandpass filter to use
bandwidth = 1;
fltsz = 31;
flt = mkBandpassCosine(res, cpim, bandwidth, fltsz, 0);
%flt = mkBandpassDog(res, cpim, bandwidth, fltsz, 0);

% Make circular stimulus masks
innerres = floor(4/totalfov * res/2)*2;
mask = makecircleimage(res,res/2*maskFrac,[],[],res/2);  % white (1) circle on black (0)

% Bars
canonicalSparsity = 5;
canonicalAngle = 0;

jumpvals = [9, 7, 5, 3, 1];
angles = [0, pi/4, pi/2, 3*pi/4];

[barsSparsity, linesSparsity, ~] = createBarStimulus(res, flt, spacing, jumpvals, nframes, canonicalAngle);

figure(2);clf; colormap gray
for ii = 1:size(barsSparsity,3)
   subplot(3,2,ii)
   imagesc(barsSparsity(:,:,ii,1)); axis image;
end
%%

histogram(output)

subplot(3,1,2:3)
imshow(output)