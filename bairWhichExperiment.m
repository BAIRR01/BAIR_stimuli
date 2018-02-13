function [experimentType, numberOfRuns] = bairWhichExperiment

% Which experiment?
experimentTypes = {'TASK' 'PRF' 'SPATIOTEMPORAL' 'HRFPATTERN'  'HRFPATTERNINVERTED'  'HRFCHECKER'  'HRFCHECKERINVERTED' };
whichExperiment = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);
experimentType = experimentTypes{whichExperiment};

% How many runs?
prompt = {'How many runs?'};
defaults = {'1'};
answer = inputdlg(prompt, 'Number of runs', 1, defaults);
numberOfRuns = str2num(answer{1,:});

end
