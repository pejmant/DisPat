% Ask the training image dimension
Dim = input('What is the training image dimension ? ') ;
Pat = input('What is the pattern dimension ? ') ;

% Read the training image
% [datain, colnames, line1]=loadgeoeas('ti2.dat');
% out = geoeas2matlab(datain,[Dim Dim]);

image=imread('ti.bmp');
image=imresize(image,[Dim Dim]);
image=im2bw(image,0.75);
out=1-image;

% show training image
% imagesc(out);

% if m=2 then every other pattern will be saved, if m=3 then 1st,4th,7th...
% will be recorded. And so on. (For reducing the dimension and redundancy).
m = input('Number of patterns to jump over in the series (Skipping) ? ');
m1 = input('Subsampling, distance between adjacent points (Multiple-Grids) ? ');
% Number of patterns : disDim*disDim
disDim = ceil((Dim - (1+(Pat-1)*m1) + 1)/m);
fprintf('\nPatterns retained for analysis = %d x %d\n\n', disDim,disDim);


X=[];
k=1;
for i=1:disDim
    for j=1:disDim
        %X((i-1)*disDim+j,:)=reshape(out(i:i+Pat-1,j:j+Pat-1),1,Pat^2);
        wx = 1+m*(i-1):m1:1+m*(i-1)+(Pat-1)*m1;
        wy = 1+m*(j-1):m1:1+m*(j-1)+(Pat-1)*m1;
       % The "if" below is to delete completely empty patterns from calculations 
       % if sum(sum(out(wx,wy)))~=0
            X(k,:)=reshape(out(wx,wy),1,Pat^2);
            k=k+1;
       % end
    end
end