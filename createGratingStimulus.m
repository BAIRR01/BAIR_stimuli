function output = createGratingStimulus(stimParams, cyclesPerDegree, gratingOrientation, peak2peakContrast)
% CREATE GRATING STIMULUS
%   stimParams - defined using stimInitialize.m; should contain a the
%       convolutional bandpass filter to use, in space domain
%       (stimParams.bpfilter)
%   stripesPerImage - Define approximate number of repetitions of lines
%       across image (previous variable name: relCutOff)
%   orientation - Orientation of grating in [0,2*pi]
%   peak2peakContrast - Michelson contrast of grating in percent [0 100];
%   Note: uses a randomly determined phase

stimWidth  = stimParams.stimulus.srcRect(3)-stimParams.stimulus.srcRect(1);
stimHeight = stimParams.stimulus.srcRect(4)-stimParams.stimulus.srcRect(2);

imSizeInDeg = stimParams.experimentSpecs.radii{1} * 2;
cpi = cyclesPerDegree * imSizeInDeg;

% Create a grating: superimpose number of orientations in
% gratingOrientation
[x,y] = meshgrid(linspace(0,1,stimWidth), linspace(0,1,stimHeight));

for ii = 1:length(gratingOrientation)
    th = gratingOrientation(ii);
    
    xr = x*cos(th) - y*sin(th);
    phX = rand*2*pi; % phase
    
    if ii == 1
        output = cos(2*pi*xr*cpi+phX) * peak2peakContrast/2/100; rms_single = rms(output(:));
    else
        output = output + cos(2*pi*xr*cpi+phX) * peak2peakContrast/2/100;
    end
end

% equalize rms contrast to that of a single grating
output = output / rms(output(:)) * rms_single;
    
return



