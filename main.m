clear all; close all; clc;
global X;global temp;global Pat;global Y;
stprpath('D:\Users\mehrdadh\stprtool');

profile on


%----------------------
bAskInputs = false;

way  = 8; 
MDS  = 10;
clus = 40;

bFastMDS = true;
%----------------------




%% Training image 
%__________________________________________________________________________

% Ask the training image dimension
Dim = input('What is the training image dimension ? ') ;
Pat = input('What is the pattern dimension ? ') ;

% Mehotd 1)Read the training image from sgems file
% [datain, colnames, line1]=loadgeoeas('ti2.dat');
% out = geoeas2matlab(datain,[Dim Dim]);

% Method 2)Read the training image from a picture
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


X=zeros(disDim*disDim,Pat^2);
k=1;
for i=1:disDim
    for j=1:disDim
        %X((i-1)*disDim+j,:)=reshape(out(i:i+Pat-1,j:j+Pat-1),1,Pat^2);
        wx = 1+m*(i-1):m1:1+m*(i-1)+(Pat-1)*m1;
        wy = 1+m*(j-1):m1:1+m*(j-1)+(Pat-1)*m1;
       % % The "if" below is to delete completely empty patterns from calculations 
       % % you should also change to X initialization to X=[];
       % if sum(sum(out(wx,wy)))~=0
            X(k,:)=reshape(out(wx,wy),1,Pat^2);
            k=k+1;
       % end
    end
end



%% Dissimilarity Matrix
%__________________________________________________________________________


disDim2=size(X,1);
if ~bFastMDS
    disMat = zeros(disDim2);
end

if bAskInputs
    disp('Which dissimilarity distance method ? ') ;
    disp(' 1 . Euclidean distance ');
    disp(' 2 . Hausdorff distance ');
    disp(' 3 . Modified Hausdorff distance ');
    disp(' 4 . Hough Transform dissimilarity distance');
    disp(' 5 . Cosine distance');
    disp(' 6 . Radon Transform dissimilarity distance');
    disp(' 7 . distance Transform ');
    disp(' 8 . distance Transform (proximity)');

    way = input(' ');
end

fprintf('\nDissimilarity Matrix number %d is being built...\n',way);

switch way
    case 8
        temp=zeros(disDim2,Pat^2);
        for i=1:disDim2
            DT=bwdist(reshape(X(i,:),Pat,Pat));
            if DT==Inf
                DT = ones(1,Pat^2);
            end
            temp(i,:)=reshape(1-DT./max(max(DT)),1,Pat^2);
        end
        if ~bFastMDS
            disMat=squareform(pdist(temp));
        end
        
    case 7
        temp=zeros(disDim2,Pat^2);
        for i=1:disDim2
            temp(i,:)=reshape(bwdist(reshape(X(i,:),Pat,Pat)),1,Pat^2);
            if temp(i,:)==Inf
                temp(i,:)= sqrt(2)*(Pat-1).*ones(1,Pat^2);
            end
        end
        if ~bFastMDS
            disMat=squareform(pdist(temp));
        end
    
    case 6
        fprintf('\n Percentage Completed:    ');
        theta=0:359;
        for i=1:disDim2
            perc=round(100*i/disDim2);
            fprintf('\b\b\b\b%3d%%', perc);
            temp(i,:)=reshape(radon(reshape(X(i,:),Pat,Pat),theta),1,[]);
        end
        if ~bFastMDS
            fprintf('\n>Calculating distances...');
            disMat=squareform(pdist(temp));
        end
        
    case 5
        fprintf('>Calculating distances...');
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            for j=i+1:disDim2
                disMat(i,j) = calcDissimilarityCosine(i,j,disDim2,Pat);
                disMat(j,i) = disMat(i,j);
            end
        end
        
     case 4
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            perc=round(100*i/disDim2);
            fprintf('\b\b\b\b%3d%%', perc);
            temp(i,:)=reshape(hough(reshape(X(i,:),Pat,Pat)),1,[]);
        end
        if ~bFastMDS
            fprintf('\n>Calculating distances...');
            disMat=squareform(pdist(temp));
        end

     case 3
        fprintf('>Calculating distances...');
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            for j=i+1:disDim2
                disMat(i,j) = calcDissimilarityModified(i,j,disDim2,Pat);
                disMat(j,i) = disMat(i,j);
            end
        end    

     case 2
        fprintf('>Calculating distances...');
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            for j=i+1:disDim2
                disMat(i,j) = calcDissimilarity(i,j,disDim2,Pat);
                disMat(j,i) = disMat(i,j);
            end
        end
        
        
    case 1     
        if ~bFastMDS
            disMat=squareform(pdist(X));
%         disMat=squareform(pdist(X),'seuclidean');
%         disMat=squareform(pdist(X),'mahalanobis');
%         disMat=squareform(pdist(X),'cityblock');
%         disMat=squareform(pdist(X),'minkowski',3);
%         disMat=squareform(pdist(X),'cosine');
%         disMat=squareform(pdist(X),'correlation');
%         disMat=squareform(pdist(X),'spearman');
%         disMat=squareform(pdist(X),'hamming');
%         disMat=squareform(pdist(X),'jaccard');
%         disMat=squareform(pdist(X),'chebychev');
        end
        
end



%% Multi Dimensional Scaling
%__________________________________________________________________________

% booleans indicating whether to plot corr. coef. or not
bPlotCorrCoefficient = true;
bPlotFinalCorrCoefficient = true;

fprintf('\n');
disp('MultiDimensional Scaling ') ;
disp('Classical  MDS');

if ~bFastMDS
    [Y1,e] = cmdscale(double(disMat));
end


% Plot EigenValues
% ---------------------------------
if ~bFastMDS
    figure;
    subplot(2,1,1);
    plot(1:min(Pat^2,disDim2),e(1:min(Pat^2,disDim2)));
    graph2d.constantline(0,'LineStyle',':','Color',[.7 .7 .7]);
    axis([1,min(Pat^2,disDim2),min(e),max(e)*0.9]);
    xlabel('Eigenvalue number');
    ylabel('Eigenvalue');
end


% Plot Correlation Coefficients
% ---------------------------------
if ~bFastMDS
    if bPlotCorrCoefficient
        subplot(2,1,2);
        hold on
        corrs=[];
        for i=1:5:min(size(Y1,2),Pat^2)
            CORR = corrcoef(disMat,squareform(pdist(Y1(:,1:i))));
            corrs=[corrs CORR(1,2)];
        end
        plot(1:5:min(size(Y1,2),Pat^2),corrs,'ko-','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerSize',8);
        xlabel('Dimension');
        ylabel('Correlation coefficient');
    end
end


% Do MDS
% ---------------------------------
if bAskInputs
    MDS = input('What is the MDS dimension ? ') ;
end

if ~bFastMDS
    maxerr = max(abs(pdist(X)-pdist(Y1(:,1:MDS)))) %#ok<NOPTS>
    Y = Y1(:,1:MDS);
else
    if size(X,1) < 100
        disMat=squareform(pdist(X));
        [Y1,e] = cmdscale(double(disMat));
        Y = Y1(:,1:MDS);
    elseif size(X,1) < 400
        [Y,totaltime] = scmdscale_withDist(x,12,50,17,1);
    else
%         [Y,totaltime] = scmdscale_withDist(X,MDS,200,floor(3*MDS/2),1);
        [Y,totaltime] = scmdscale_withDistAndEigs(X,MDS,200,floor(3*MDS/2),1);
    end
end

% Plot final distances' correlation
% ---------------------------------
if ~bFastMDS
    if bPlotFinalCorrCoefficient
        MDStest = squareform(pdist(Y));
        figure;
        plot(disMat,MDStest,'b.');
        CORR = corrcoef(disMat,MDStest);
        xlabel('Dissimilarity distance');
        ylabel(['Euclidean distance in ',num2str(MDS),' D']);
        title(['Correlation Coefficient = ',num2str(CORR(1,2))]);
    end
end



%% clustering Algorithms
%__________________________________________________________________________


bShow3DProjection = true;

% Choosing between clustering methods
disp('Clustering Method ? ') ;
disp('Kernel K-means analysis ');
way2 = 3;
MDS1=MDS;

if bAskInputs
    clus = input('How many clusters ? ') ;
end

sigma =0.2*max(max(pdist(Y)));

% Apply Radial Basis Function Kernel
K = kernel(Y', 'rbf', sigma);
% % Faster Kernel about 15%
% K = slkernel(Y', 'gauss', sigma);


% Spectral Clustering to find initial clusters
idx = SpectralClustNg( K, clus);


% Define the dimension of the projections on the Feature space
% e = eig(K);
%new_dim=length(find(e>1));
new_dim=3;
model= kpca(Y', K, struct('ker','rbf','arg',sigma,'new_dim',new_dim));
Z = kernelproj(Y', model)';


% Kernel K-Means algorithm with initial points from Spectral Clustering.
idx = dualkmeansFast(K,clus,idx);


MDS=-1;




%% plotting the results
%__________________________________________________________________________

figure;
ptsymb = {'bs','r^','md','go','c+','y*','kv','b>','r<','mp','gh','c.'};
ptsymb=repmat(ptsymb,1,8); PC=0;


for i = 1:clus
    clust = find(idx==i);
    plot3(Z(clust,1),Z(clust,2),Z(clust,3),ptsymb{i});
    hold on
end

title('KPCA results in 3-D kernel space with centers');
grid on;


%% showing average of all clusters patterns
%__________________________________________________________________________


bShowAverageOfClusters = false;

if bShowAverageOfClusters
    calcClusters
end

MDS = MDS1;



%% Plotting the clustered categories on training image wrt. their colors
%__________________________________________________________________________


figure;
Prototype=zeros(Dim,Dim);
for i = 1:clus
    clust = find(idx==i);
    Iclust = ceil(clust./disDim);
    Jclust = mod(clust,disDim)+1;
    for j = 1:length(Iclust)
        Prototype(1+m.*(Iclust(j)-1)+(Pat-1)*m1/2,1+m.*(Jclust(j)-1)+(Pat-1)*m1/2)=i;
    end
end
imagesc(Prototype);
axis square



%% Find Random Point in a Cluster

% Initialization of the model
clusterModel = kpca(Y', K, struct('ker','rbf','arg',sigma,'new_dim',Pat^2));
% clusterModel = greedykpca(Y',struct('ker','rbf','arg',sigma,'new_dim',Pat^2)); 
Z = kernelproj(Y', clusterModel)';

% Finding a Point in cluster using the model
cluster = clusterStatistics(Z, 10, idx);
randPattern = lhs_iman_n(cluster.mean, cluster.std, cluster.corr);

randPattern_inR = Yrec(randPattern', clusterModel)';



%% Pattern warping
% using Radon Algorithm
[distance,index] = find_closest_pattern_index(randPattern_inR,3);
[I,threshold] = Radon_warp(distance,index);
Pattern = I>threshold;
figure;drawPattern(Pattern);



profile off
profile viewer
