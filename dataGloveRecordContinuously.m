[subjID, sessionID, ssDefined] = bairWhichSubjectandSession();
if ~ssDefined, return; end

glovePointer = initializeDataGlove();

t0 = GetSecs();
sampleTime = 2/60;
counter = 1;

pth = fullfile(vistadispRootPath, 'Data');
fname = sprintf('sub-%s_ses-nyumeg%2d_DataGlove_%s',subjID, sessionID, datestr(now, 'yyyy_mm_dd_hh_MM_ss'));

fid = fopen(fullfile(pth, fname) ,'w');

f = figure; hold on
for ff = 1:5
    f(ff) = animatedline;
end
recordData = true;

while recordData
    recordTime = sampleTime * counter;
    t = GetSecs() - t0;
    WaitSecs(recordTime - t);
    data = sampleDataglove (glovePointer);
    
    
    fprintf(fid,  '%s\t', datestr(now, 'hh_MM_ss_FFF'));
    for ii = 1:length(data)
        fprintf(fid,  '%d\t', data(ii));
    end
    fprintf(fid, '\n');
    for ff = 1:5
        addpoints(f(ff), counter, data(1,ff))
    end
    drawnow
    counter = counter + 1;
    
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
end

fclose(fid);
