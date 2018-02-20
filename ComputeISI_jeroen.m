% WE ONLY USE STIMULUS TRAIN OF RUN1, other runs are not used

D=load('NYU_run1.mat');
TIMING_NYU=D.VarName1';
ISI_NYU=diff(TIMING_NYU);
TR=850;
factor=5;
TReff=TR/factor;
jitteroffset_ms=TReff;
% ISI in TRs, useing NYU ISI distribution as a start
ISI_jeroen=round(ISI_NYU*1000./TR);
%ISI_jeroen=[6,20,13,13,22,5,10,6,8,6,8,15,5,9,16,18,4,14,8,5,11,9,7,17,29,5,10,12,7,4,25];

% use randi(5,[1,32])-ones(1,32) to make jitter offsets equal to 1/5 of TR,
% in this case 850/5 = 170ms, our effective TR
ISI_jeroen_jittoffsets=[0,1,2,1,3,4,1,3,0,4,0,2,1,4,2,2,4,1,2,4,0,4,0,3,2,3,1,3,4,3,0];
%ISI_jeroen_jittoffsets=zeros(1,31);
negpos_Array=[-1,1,1,1,-1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,1,-1,1,1,1,1];
%some jitter offsets should be subtracted
ISI_jeroen_jittoffsets=ISI_jeroen_jittoffsets.*negpos_Array;

%now convert everything to ms, where the idea is that the stimulus will
%always coincide with a slice acquisition (will also work for multiband)

ISIjitter_jeroen_ms=ISI_jeroen*TR + ISI_jeroen_jittoffsets*jitteroffset_ms;
ISIjitter_jeroen_s=ISIjitter_jeroen_ms/1000;

%convert back to timing array, just as TIMING_NYU
%  BASELINE 
baseline=36; %integer number of TRs
TIMING_jeroen_s=cumsum([baseline*TR/1000,ISIjitter_jeroen_s])';
TIMING_jeroen_ms=TIMING_jeroen_s*1000;

%save stimuli UMCU, NYU
FID=fopen('ISI_UMCU_ms.txt','w');
fprintf(FID, '%6.0f \n',TIMING_jeroen_ms);
fclose(FID);
FID=fopen('ISI_NYU_s.txt','w');
fprintf(FID, '%3.3f \n',TIMING_jeroen_s);
fclose(FID);

figure, 
subplot 121,hist(ISI_NYU,10), title('ISI NYU ORIGINAL')
subplot 122,hist(ISIjitter_jeroen_s,10), title('ISI modified to match slicetiming')

disp(num2str(TIMING_jeroen_s))