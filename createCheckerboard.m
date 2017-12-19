function im = createCheckerboard(stimParams, cpi)
% CREATE checkerboard stimulus
%   stimParams - structure created by stimInitialize
%   cpi    - size of checks (cycles per image)

stimWidth  = stimParams.stimulus.srcRect(3)-stimParams.stimulus.srcRect(1);
stimHeight = stimParams.stimulus.srcRect(4)-stimParams.stimulus.srcRect(2);

[x,y] = meshgrid(linspace(0,1,stimWidth), linspace(0,1,stimHeight));


th = pi/4;

xr = x*cos(th) - y*sin(th);
yr = x*sin(th) + y*cos(th);

phX = rand*2*pi;
phY = rand*2*pi;
plaid = cos(2*pi*xr*cpi+phX)+cos(2*pi*yr*cpi+phY);

im = double(plaid>0)-.5;


% figure, 
% imshow(im+.5)

return

end

