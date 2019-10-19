cd('/Users/winawerlab/Documents/temp/sixcatlocalizer');
load('sixcatlocalizer')
%%
for ii = 171:size(bodies,4)
     I = bodies(:,:,:,ii);
    Imode = mode(I(:));
    if Imode ~= 128
       ind = sum(I==Imode,3)>1;
       %I(I == Imode) = backgroundColor;
       for dim = 1:size(I,3)
           temp_I = I(:,:,dim);
           temp_I(ind) = 128;
           I(:,:,dim) = temp_I;
       end
    end   
    figure;imshow(I);
    title(ii)
    set(gcf, 'Position', [595         304        1195        1041] ); 
    waitforbuttonpress;close;
end