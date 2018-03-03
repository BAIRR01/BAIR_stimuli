function experimentType = bairWhichExperimentList

% Which experiment?
whichExperiment = listdlg('PromptString', 'Which experiment?', 'ListString', experimentTypes);
experimentType = experimentTypes{whichExperiment};

end
