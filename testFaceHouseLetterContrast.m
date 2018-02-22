load('/Volumes/server/Projects/BAIR/Stimuli/Kay2013_unfilteredstimuli/stim.mat');
load('/Users/winawerlab/matlab/toolboxes/vistadispBAIR/StimFiles/spatiotemporal_NYU-ECOG_102.mat')

% find images 
% 
% %%
% for ii = 1:size(houses,3)
%     figure;imshow(houses(:,:,ii),[])
%     title(ii);
%     waitforbuttonpress;close;
% end

%faces
figure;hold on;
pli = 1;
for ii = 73:80%96%1:size(stimulus.images,3)
    subplot(2,8,pli);
    imshow(stimulus.images(:,:,ii))
    title(ii);
    pli = pli+1;axis tight;
end

Flist = [30 35 24 83 74 13 59 41];
for ii = 1:8%96%1:size(stimulus.images,3)
    subplot(2,8,pli);
    imshow(faces(:,:,Flist(ii)))
    title(Flist(ii));
    pli = pli+1;axis tight;
end

% letters
figure;hold on
pli = 1;
for ii = 81:88%96%1:size(stimulus.images,3)
    subplot(2,8,pli);
    imshow(stimulus.images(:,:,ii))
    title(ii);
    pli = pli+1;axis tight;
end

Llist = [13 21 6 20 1 9 8 24];
for ii = 1:8%96%1:size(stimulus.images,3)
    subplot(2,8,pli);
    imshow(letters(:,:,Llist(ii)),[])
    title(Llist(ii));
    pli = pli+1;axis tight;
end

% houses
figure;hold on
pli = 1;
for ii = 89:96%96%1:size(stimulus.images,3)
    subplot(2,8,pli);
    imshow(stimulus.images(:,:,ii))
    title(ii);
    pli = pli+1;axis tight;
end

Hlist = [33 7 24 20 41 3 31 25];
for ii = 1:8%96%1:size(stimulus.images,3)
    subplot(2,8,pli);
    imshow(houses(:,:,Hlist(ii)),[])
    title(Hlist(ii));
    pli = pli+1;axis tight;
end

%%
load('/Volumes/server/Projects/BAIR/Stimuli/Kay2013_unfilteredstimuli/stim.mat');
I0 = faces(:,:,30); %houses(:,:,33); 

% resize such that FHL are ~25% of pattern stimuli 
I1 = imresize(I0, 0.4*imageSizeInPixels); 
stimSizeInPixels = size(I1);

% paste into a grayscale background of imageSizeInPixels
I2 = padarray(I1,[(imageSizeInPixels(1)-size(I1,1))/2 (imageSizeInPixels(1)-size(I1,1))/2], mode(I1(1,:)));

% band-pass filter
I3 = imfilter(I2, stimParams.bpFilter);

% apply circular mask to get rid of edges and background noise
%supportDiameter       = imageSizeInPixels;
%maskRadius            = 1.1*(1*size(I1,1))/2;
%maskOrigin            = [(imageSizeInPixels(1)+1)./2+(0.05*imageSizeInPixels(1)) (imageSizeInPixels(1)+1)./2 ];
%maskOrigin            = [(imageSizeInPixels(1)+1)./2 (imageSizeInPixels(1)+1)./2];
%circularMask          = mkDisc(supportDiameter, maskRadius, maskOrigin, 1/100 * imageSizeInPixels(1), [1 0]);
maskRadius             = 0.8*stimSizeInPixels(1)/2;
maskOrigin             = [(imageSizeInPixels(1)+1)./2 (imageSizeInPixels(1)+1)./2];

squareMask             = makecircleimage(imageSizeInPixels(1),maskRadius,[],[], 1.1*maskRadius,1) .* ...
                            makecircleimage(imageSizeInPixels(1),maskRadius,[],[], 1.1*maskRadius,2);
maskRadius             = 0.8*stimSizeInPixels(1)/2;

circularMask           = makecircleimage(imageSizeInPixels(1),maskRadius,[],[],1.1*maskRadius,0, [maskOrigin(1)+0.1*maskOrigin(1) maskOrigin(2)], [1 1.2]);
combinedMask           = circularMask .* squareMask;
%I4                     = bsxfun(@times,I3,squareMask);
%I5                     = bsxfun(@times,I3,circularMask);
I4                     = bsxfun(@times,I3,combinedMask);

%figure;imshow(I4)
%figure;imshow(I5)
%figure;imshow(I6);
%figure;imshow(I4-I5);

% scale contrast to prevent (too much) clipping
%scaleFactor = range(prctile(I4(I4~=0), [2 98]));
scaleFactor = 1.5;
I5 = I4 * scaleFactor/range(I4(I4~=0));

% 8Bit
I6 = uint8((I5+.5)*255);


%figure;
subplot(3,4,1);imshow(I0);axis on; title('original');
subplot(3,4,2);histogram(I0);
set(gca, 'XLim', [0 1], 'YLim', [0 3*10^4]);

subplot(3,4,3);imshow(I2);axis on; title('resize & fit in master display size');
subplot(3,4,4);histogram(I2);
set(gca, 'XLim', [0 1], 'YLim', [0 2*10^4]);

subplot(3,4,5);imshow(I3);axis on; title('bp filter');
subplot(3,4,6);histogram(I3);
set(gca, 'XLim', [-5 5], 'YLim', [0 2*10^4]);

subplot(3,4,7);imshow(I4);axis on; title('circular mask');
subplot(3,4,8);histogram(I4);
set(gca, 'XLim', [-5 5], 'YLim', [0 2*10^4]);

subplot(3,4,9);imshow(I5);axis on;title('scale (to prevent too much clipping)');
subplot(3,4,10);histogram(I5(I5~=mode(I5(:)))); title('(mode excluded from hist)');
set(gca, 'XLim', [-1 1])% 'YLim', [0 0.1*10^4]);

subplot(3,4,11);imshow(I6);axis on;title('8bit');
subplot(3,4,12);histogram(I6(I6~=mode(I6(:)))); title('(mode excluded from hist)');
set(gca, 'XLim', [0 255]); %'YLim',[0 0.1*10^4]);

         
subplot(3,4,3);imshow(I2);axis on; title('resize & fit in master display size');
subplot(3,4,4);histogram(I2);
set(gca, 'XLim', [0 1], 'YLim', [0 2*10^4]);

subplot(3,4,5);imshow(I3);axis on; title('bp filter');
subplot(3,4,6);histogram(I3);
set(gca, 'XLim', [-5 5], 'YLim', [0 2*10^4]);

subplot(3,4,7);imshow(I4);axis on; title('circular mask');
subplot(3,4,8);histogram(I4);
set(gca, 'XLim', [-5 5], 'YLim', [0 2*10^4]);

subplot(3,4,9);imshow(I5);axis on;title('scale (to prevent too much clipping)');
subplot(3,4,10);histogram(I5(I5~=mode(I5(:)))); title('(mode excluded from hist)');
set(gca, 'XLim', [-1 1])% 'YLim', [0 0.1*10^4]);

subplot(3,4,11);imshow(I6);axis on;title('8bit');
subplot(3,4,12);histogram(I6(I6~=mode(I6(:)))); title('(mode excluded from hist)');
set(gca, 'XLim', [0 255]); %'YLim',[0 0.1*10^4]);

         