function rootPath=BAIRRootPath()
% Return the path to the root BAIR stimulus directory
%
% This function must reside in the directory at the base of the BAIR
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(BAIRRootPath,'stimuli')

rootPath=which('BAIRRootPath');

rootPath=fileparts(rootPath);

return
