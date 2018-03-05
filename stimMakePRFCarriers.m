function barCarriers = stimMakePRFCarriers(stimParams)

imageSizeInPixels = size(stimParams.stimulus.images);

site = stimParams.experimentSpecs.Row{1};

files = dir(fullfile(vistadispRootPath, 'StimFiles', ...
    sprintf('%s_hrfpattern*', site)));

for ii = 1:length(files)
    load(fullfile(files(ii).folder, files(ii).name), 'stimulus');
    n = size(stimulus.images,3);
    im = stimulus.images(:,:,1:(n-1)/2);
    if ii == 1, images = im; else, images = cat(3, images, im); end
end

barCarriers = imresize(images, imageSizeInPixels, 'nearest');
barCarriers = cat(3, barCarriers, flip(barCarriers,1));
barCarriers = cat(3, barCarriers, flip(barCarriers,2));

return