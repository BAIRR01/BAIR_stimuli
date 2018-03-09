function [numberOfExperiments, experimentTypes, runIDs, selectionMade] = bairWhichExperimentList(site)

% Which experiment list file?
filespec = fullfile(vistadispRootPath, 'RunMe', '*.txt');
[fname, pathname, fileselected] = uigetfile(filespec, 'Select a text file with the experiment list');

if fileselected
    % Read the content of the selected text file
    fileID = fopen(fullfile(pathname,fname));
    TXT = textscan(fileID, '%s');
    fclose(fileID);

    % Parse the filename string
    TXT = TXT{:};
    numberOfExperiments = size(TXT,1);

    for ii = 1:numberOfExperiments
        thisExperimentName = TXT{ii,:};
        dashPos = strfind(thisExperimentName,'_');
        experimentTypes{ii} = thisExperimentName(1:dashPos-1);
        runIDs(ii) = str2num(thisExperimentName(dashPos+1:end));
    end

    % DEBUG: code below doesn't work when file has just a single experiment
    % T = readtable(fullfile(pathname, fname));
    % experimentType      = T.Var1;
    % runID               = T.Var2;
    % numberOfExperiments = height(T);

    % Check that the experiment files exist
    stimPath = fullfile(vistadispRootPath, 'StimFiles');
    for ii = 1:numberOfExperiments
        fname = sprintf('%s_%s_%d.mat', site, experimentTypes{ii}, runIDs(ii));

        if ~exist(fullfile(stimPath, fname), 'file')
            error('Requested experiment file %s not found in expected location:\n%s', fname, stimPath);        
        end
    end

    fprintf('All experiment files ready to go!\n');
    selectionMade = 1;
else
    numberOfExperiments = [];
    experimentTypes = [];
    runIDs = [];
    selectionMade = 0;
end

end
