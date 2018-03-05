experiments = { ...
    'hrfpattern_1' ...
    }';

experimentList = table(experiments);

fname = 'experimentsToRun';
pth = vistadispRootPath;
writetable(experimentList, fullfile(pth, sprintf('%s.txt', fname)), ...
            'FileType','text', 'Delimiter', '\t', 'WriteVariableNames', false)
        
        
%fileId = fopen('experimentsToRun.txt', 'w')
%fprintf