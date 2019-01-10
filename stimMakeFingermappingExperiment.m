function stimMakeFingermappingExperiment(stimParams,  runNum, TR, experimentType)

path = fullfile(sensorimotorRootPath, 'motor/UMCU-Code/FingerMapIEMU/01_SeqRightHand');

%load timing from C++ code directory
timing = load(fullfile(path, 'timing.dat'));

%finger flex pattern from C++ code directory
blockInfo = load(fullfile(path, 'blockinfo.dat'),'-ascii');