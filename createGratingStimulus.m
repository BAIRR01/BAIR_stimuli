function output = createGratingStimulus(stimParams, cyclesPerDegree, scaleContrast, stripeAngle)
% CREATE GRATING STIMULUS
%   stimParams - defined using stimInitialize.m; should contain a the
%       convolutional bandpass filter to use, in space domain
%       (stimParams.bpfilter)
%   stripesPerImage - Define approximate number of repetitions of lines
%       across image (previous variable name: relCutOff)
%   scaleContrast - Scale the edge contrast relative to the
%       empirically found maximum that prevented too much clipping for SOC
%       stimuli (default = maximum);
%   stripeAngle - Orientation of grating in [0,2*pi]
%   Note: uses a randomly determined phase

if isempty(scaleContrast)
    % Scale the edge contrast (this value was found empirically to maximize
    % contrast whilst preventing too much clipping for SOC stimuli)
    scaleContrast = 0.18;
else
    scaleContrast = scaleContrast * 0.18;
end

stimWidth  = stimParams.stimulus.srcRect(3)-stimParams.stimulus.srcRect(1);
stimHeight = stimParams.stimulus.srcRect(4)-stimParams.stimulus.srcRect(2);

imSizeInDeg = stimParams.experimentSpecs.radii{1} * 2;
cpi = cyclesPerDegree * imSizeInDeg;

% Create a grating
[x,y] = meshgrid(linspace(0,1,stimWidth), linspace(0,1,stimHeight));

for ii = 1:length(stripeAngle)
    th = stripeAngle(ii);
    
    xr = x*cos(th) - y*sin(th);
    phX = rand*2*pi; % phase
    
    if ii == 1
        output = cos(2*pi*xr*cpi+phX) / 2; rms_single = rms(output(:));
    else
        output = output + cos(2*pi*xr*cpi+phX) / 2;
    end
end

output = output / rms(output(:)) * rms_single;
return



