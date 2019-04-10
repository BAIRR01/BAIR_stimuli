function generate_stimuli(condition, scan_period, isi_min, isi_max, max_dur)
%GENERATE_STIMULI write stimuli (onsets and bitmaps) in the correct scheme 
% folder
%
% INPUTS:
%   condition : str
%     one of 'practice', 'IEMU', 'fMRI' (case is important)
%   scan_period : int
%     should be identical to the value of "scan_period" in the Presentation
%     file
%   isi_min : int
%     mininum inter-stimulus distance, between onset times, in s
%   isi_max : int
%     maximum inter-stimulus distance, between onset times, in s
%   max_dur : int
%     maximum duration of the whole task, in s
%
% EXAMPLES:
%   generate_stimuli('practice', 1000, 4, 6, 300)
%
%   generate_stimuli('fMRI', 850, 6, 15, 480)
%
%   generate_stimuli('IEMU', 1000, 4, 6, 480)
%
% gio@gpiantoni.com  2018/09/26

n_bitmaps = 4;
port_code_offset = 10; 

% convert from seconds to scans
isi_min1 = round(isi_min / scan_period* 1000);
isi_max1 = round(isi_max / scan_period* 1000);
max_dur1 = round(max_dur / scan_period* 1000);
fprintf('ISI (in integer scans), min: %d, max: %d, max number of scans: %d\n', ...
    isi_min, isi_max, max_dur);

% PATHS 
dir_schemes = fileparts(mfilename('fullpath'));
dir_name = ['scheme_gestures_' condition];
dir_output = fullfile(dir_schemes, dir_name);
 
file_onsets = fullfile(dir_output, 'picture_onset_sequence.txt');
file_bitmaps = fullfile(dir_output, 'bitmap_filename_sequence.txt');
file_codes = fullfile(dir_output, 'picture_port_code_sequence.txt');

% NUMBER OF EVENTS
n_events = floor(max_dur / ((isi_max + isi_min) / 2));
fprintf('Number of events: %d\n', n_events)

% ONSETS
isi = [1; randi(isi_max - isi_min + 1, n_events - 1, 1) + isi_min - 1];
onsets = cumsum(isi);

s = sprintf('%d\n\r', onsets);
fid = fopen(file_onsets, 'wt');
fprintf(fid, s(1:end-2));
fclose(fid);

% EVENTS
events = randi(n_bitmaps, n_events, 1);

fid = fopen(file_bitmaps, 'wt');
s = sprintf('exec_stim_%d.jpg\n\r', events);
fprintf(fid, s(1:end-2));
fclose(fid);

s = sprintf('%d\n\r', events + port_code_offset);
fid = fopen(file_codes, 'wt');
fprintf(fid, s(1:end-2));
fclose(fid);