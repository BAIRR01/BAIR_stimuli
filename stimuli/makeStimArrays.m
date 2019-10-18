clear all;
%stimDir = '/Users/iris/fromNIHcomputer/groenii/Experiments_local/6categorylocalizer/stims';
stimDir = '/Users/iris/DataAnalyses_local/IRIS_H6_WeibullStop/tmp_stimuli';

cd(stimDir);
d = dir(stimDir);

for ii = 1:length(d)
    if d(ii).isdir && ~contains(d(ii).name, '.')
        varname = genvarname(d(ii).name);
        disp(varname);
        catDir = fullfile(d(ii).folder, d(ii).name);
        d2 = dir(fullfile(catDir, '*.jpg'));
        I = imread(fullfile(catDir, d2(1).name));
        stimArray = nan([size(I) length(d2)]);
        for jj = 1:length(d2)
            stimArray(:,:,:,jj) = imread(fullfile(catDir, d2(jj).name));
        end
        stimArray = uint8(stimArray);   
        eval([varname '= stimArray;' ]);            
    end        
end

clear I d d2 ii jj stimDir catDir varname stimArray
save('sixcatlocalizer');

