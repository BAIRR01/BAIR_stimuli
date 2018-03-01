function runID = bairWhichRun

% Which run
prompt = {'Which run?'};
defaults = {'1'};
answer = inputdlg(prompt, 'Run ID', 1, defaults);
runID = str2num(answer{1,:});

end
