% Records from the dataglove continuously (for use at MEG)

% get the subject and session for file name
[subjID, sessionID, ssDefined] = bairWhichSubjectandSession();
if ~ssDefined, return; end

% initialize the dataglove and set t0
glovePointer = initializeDataGlove();
t0           = GetSecs();
sampleTime   = 2/60;
counter      = 1;
recordData   = true;

% set path and file name for saving
pth       = fullfile(vistadispRootPath, 'Data');
printTime = datestr(now, 'yyyy_mm_dd_hh_MM_ss');
fname     = sprintf('sub-%s_ses-nyumeg%s_DataGlove_%s',subjID, sessionID, printTime);
fid       = fopen(fullfile(pth, fname) ,'w');

% initialize figure and a line for each finger
f = figure; hold on
for ff = 1:5
    f(ff) = animatedline;
end

% Start recording
while recordData
    recordTime = sampleTime * counter;
    t          = GetSecs() - t0;
    % wait until specific times to record data 
    WaitSecs(recordTime - t);
    data = sampleDataglove (glovePointer);
    
    % write it to the file
    fprintf(fid,  '%s\t', datestr(now, 'hh_MM_ss_FFF'));
    for ii = 1:length(data)
        fprintf(fid,  '%d\t', data(ii));
    end
    fprintf(fid, '\n');
    
    % plot it 
    for ff = 1:5
        addpoints(f(ff), counter, data(1,ff))
    end
    drawnow
    
    % check for a q as keyboard input to stop recording
    [ssKeyIsDown,~,ssKeyCode] = KbCheck(-1);
    
    if ssKeyIsDown
        str = KbName(find(ssKeyCode));
        if iscell(str)
            str = str{1};
        end
        str = str(1);
        if str == 'q'
            recordData = false;
        end
    end
    counter = counter + 1;
end

fclose(fid);
