glovePointer = initializeDataGlove();

t0 = GetSecs();  
sampleTime = 2/60;
counter = 1;

pth = fullfile(vistadispRootPath, 'Data');
fname = sprintf('subj00DataGlove_%s', datestr(now, 'yyyy_mm_dd_hh_MM_ss'));

fid = fopen(fullfile(pth, fname) ,'w'); 

while true
    recordTime = sampleTime * counter;
    t = GetSecs() - t0;
    WaitSecs(recordTime - t);

    data = sampleDataglove (glovePointer);
    
    fprintf(fid,  '%s\t', datestr(now, 'hh_MM_ss_FFF'));
    for ii = 1:length(data)
       fprintf(fid,  '%d\t', data(ii));      
    end
    
    fprintf(fid, '\n');
 
    counter = counter + 1;
end

fclose(fid);
