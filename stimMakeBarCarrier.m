function bar_carrier = stimMakeBarCarrier(stimParams)


[output, edge, thresh, res] = createPatternStimulus([768, 768], 1/40, stimParams.bpfilter);

% Download bar stimuli
bar_carrier = [];
for ii = 1:16
    readPth = sprintf('https://wikis.nyu.edu/download/attachments/85394548/bar_carrier%d.mat?api=v2', ii);
    stimDir = fullfile(BAIRRootPath, 'stimuli');
    fname = sprintf('bar_carrier%d.mat', ii);
    writePth = fullfile(stimDir, fname);
    if ~exist(writePth, 'file'), websave(writePth,readPth); end
    im = load(writePth);
    f = fieldnames(im);
    bar_carrier = cat(3, bar_carrier, im.(f{1}));
end
