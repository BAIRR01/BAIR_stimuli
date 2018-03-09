 
load  Master_hrfpattern_1.mat
 
%load  UMC-7T_spatialpattern_1.mat
 
% onsets in ms
tr = 850; % ms
%onsets = stimulus.tsv.onset*1000;
onsets = stimulus.seqtiming(stimulus.seq~=65)*1000;
n = length(onsets);
 
figure, 
stem(mod(onsets,tr), ones(1,n))

xlim([0 850])


figure;
subplot(2,1,1);
stem(stimulus.seqtiming(stimulus.seq~=65), ones(length(stimulus.seqtiming(stimulus.seq~=65)),1), 'b')
subplot(2,1,2);
stem(stimulus.tsv.onset, ones(length(stimulus.tsv.onset),1), 'b')

figure;
subplot(2,1,1);
stem(stimulus.tsv.onset, stimulus.tsv.stim_file_index, 'b')
subplot(2,1,2);
stem(stimulus.seqtiming(stimulus.seq~=65), stimulus.seq(stimulus.seq~=65), 'b')

load  UMC-7T_spatialpattern_1.mat

figure;
subplot(3,1,1);
stem(stimulus.tsv.onset, stimulus.tsv.stim_file_index, 'b')
xlim([0 200])
subplot(3,1,2);
stem(stimulus.seqtiming_sparse(stimulus.seq_sparse~=37), stimulus.seq_sparse(stimulus.seq_sparse~=37), 'b')
xlim([0 200])
subplot(3,1,3);
stem(stimulus.seqtiming(stimulus.seq~=37), stimulus.seq(stimulus.seq~=37), 'b')
xlim([0 200])

load('/Users/winawerlab/Box Sync/Stimuli/FOR UMC/UMC-7T_temporalpattern_1.mat')

max(rem(stimulus.onsets,0.170))% max distance of onsets from 1/5 of TR
