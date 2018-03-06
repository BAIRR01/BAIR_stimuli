function experimentType = bairWhichExperiment

% Which experiment?
experimentTypes = {'TASK' 'PRF' 'SPATIALPATTERN' 'SPATIALOBJECT' ...
    'TEMPORALPATTERN' 'HRFPATTERN' 'HRFPATTERNBREATHINGCHALLENGE'};
% 'HRFCHECKER'  'HRFPATTERNINVERTED' 'HRFCHECKERINVERTED' };
whichExperiment = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);
experimentType = experimentTypes{whichExperiment};

end
