DisPat
======

MATLAB codes for DisPat
Please do the following to use the codes:
1.  Copy-paste the DisPat codes in a specific directory of your computer.
2.	Go to ...\DisPat\dist_pat folder and open distPat.m file.
3.	Press f5 or run to execute the codes.
4.	If the matlab editor asked to change the folder or add to path, please select change folder. It adds all of the necessary paths to your matlab path.
5.	You can change different parameters such:

par.Dim  : Size of training image in X and Y directions;
par.Dimz : Size of training image in Z direction;
par.Pat  : Size of template in X and Y directions;
par.Patz : Size of template in Z direction;
par.innerPatch      : Size of inner patch in X and Y directions;
par.innerPatchz     : Size of inner patch in Z direction;
par.multipleGrid    : Number of multi-grid;
 
par.bShowMultiGrids : Flag to whether show the multi-grids results or not;
 
par.hardData : Flag to whether it is conditional simulation or not;
 
 
% histogram transformations
par.bTransCat    : Flag to whether transfer the histogram of categorical training images or not;
par.bTransCon    : Flag to whether transfer the histogram of continuous training images or not;
 
 
par.m    : Number of patterns to be skipped;
par.way  : The distance function for similarity measurement;
par.MDS  : Number of dimension of MDS space;
par.clus : 100; Number of clusters
par.bUseKernelForClustering  : Flag to whether use kernel k-means or not;
par.bSkipPreviousFrozenNodes : Flag to ask about skipping from previous frozen nodes;
par.bUseDualTemplate         : Flag to ask whether to use dual template or not;
 
 
% Load/Save files
par.bLoadVariables : Load related parameters from saved files;
par.bSaveVariables : Save related parameters into files;
 


getTrainingImage 
Gets the dimension of the training image and pattern template, with skip-size and multi-grids, and then reads the file and stores the patterns in X (pattern database) each row represents a pattern. Note, in order to see Prototype, Dim and Pat should be odd values.

distanceMatrix 
Calculate the dissimilarity matrix with different methods using X. It puts the result in disMat matrix.

multiDimScale 
Runs multidimensional scaling with 4 different methods. The main one is the first method, where it first does cmdscale and then plots the eigenvalues and correlation coefficients , and then asks for the MDS dimension, and puts the results in Y.

patCluster 
Pattern clustering, with 3 methods, where the first and third one are important. In the third one, it used kernel k-means and puts the number of clusters in clus and uses dualkmeans3 file to do the clustering. It also stores the kernel in 3d in Z.

plotPoints 
Plots the results with different clusters in different colors. For the case of kernel where MDS is -1 it plots also the mesh.

showPatterns 
If the dimensionality of MDS was 2, by clicking on the plot of points you will get the index of the corresponding pattern.

showClusters 
Shows the average of all cluster patterns. It uses calcClusters to do it.

showPrototype 
Plots the clustered categories on the training image.

Functions


File and TI manipulations

matlab2geoeas(DataIn)
Transform Matlab array into GeoEAS column vector

loadgeoeas(file)
load geoeas formatted file, GSLIB formatted file.


geoeas2matlab(datain,[Dim Dim])
Reads geos to Matlab variable, it first has to be loaded with loadgeoeas function. Like example below:
[datain, colnames, line1]=loadgeoeas('ti2.dat');
out = geoeas2matlab(datain,[Dim Dim]);


Dissimilarity Matrix Calculations

calcDissimilarity(i,j,disDim2,Pat)
calculates the Hausdorff distance between X(i) and X(j) patterns, disDim2 is only needed for progressbar. And Pat is the size of the template.

calcDissimilarityCosine(i,j,disDim2,Pat)
calculates the cosine dissimilarity matrix (not working very well)

calcDissimilarityModified(i,j,disDim2,Pat)
calculates the Modified Hausdorff distance between X(i) and X(j) patterns

distanceTest
It uses patNum=[1258 123 234]; where pattern indexes are given and then calculates the distances between them and plots them. It is a test file for comparison purposes only.


MDS methods

IsomapII

L2_distance(a,b,df)
Computes Euclidean distance matrix

lle(X,K,d)
LLE algorithm for dimensionality reduction

valid_findk
valid_index_plot
valid_internal deviation
valid_sumsquares
valid_index
ind2cluster(labels)
xlim(arg1, arg2)
Used for validity checking on the number of clusters in kmeans algorithm. Where validity_Index is used to do everything.


Clustering methods

kmeans2(X, k, varargin)
K means algorithm

Dualkmeans3
K-means in kernel space

distpoint_centroid(ff,K,N)
Function calculating the distance between the centroid and the points of the cluster

MeanShiftCluster(dataPts,bandWidth,plotFlag)
Mean shift clustering method (not used anymore)


Kernel analysis

featuredist(i,j,K)
Compute the distance between two points in the feature space. Using kernel K.

SpectralClustNg(Y,sigma,N)
Peforms spectral K-means for all samples specified by the kernel K. It uses kmeans2 and centroidortho2 functions.

centroidortho2

Y = Yrec(Z,model,R_dim)
Z is the projection in kernel space, and it reconstructs the points in R input pace.


Cluster analysis

KMDL(idx,K,d)
Finds KMDL where idx holds class labels and K is kernel and d is the size of MDS space (size(Y,2))

findKMDL
calculates the KMDL and plots the results. Maximum values corresponds to the best number of clusters. It uses KMDL function to calculate KMDL value at each Nb of cluster.

calcClusters
calculate the prototypes and shows them. In here you need X, clus (number of clusters) and idx (holding the corresponding cluster that each Xi belongs to). It calculates four different ways of sharpness indexes.

sgemsclusters
finds sgems clusters prototype plot, for the case of 51 Dim and 9 Pat. Using prot51.dat where ti.bmp is used as the reference one.

showPatternsInCluster(i,idx)
Shows the patterns in cluster i according to class labels in idx.

showPrototype
shows the prototype of the TI clustrered set. And plots it.

clustercenter
calculates the center of a cluster in the kernel projected space.


Random Sampling help

clusterStatistics(i,idx)
I being the cluster number and idx being the class labels matrix. We also need Y, sigma, Pat which are global here. We model it to Pat2 dimensions and kernel project it, to find cluster.mean, cluster.std, cluster.corr, cluster.model.

latin_hs(xmean,xsd,nsample,nvar)
latin hypercube sampling with no correlation

LHS
Uses LHS on the circle data and generates the plots shown in affiliate paper of 2008.

lhs_iman_n(xmean,xsd,corr)
LHS with correlation, this is a good one as it uses mchol function to find the cholesky decomposition.

ransamp(xmean,xsd,corr,nsample)
Random sampling with correlation (not good enough because of chol function non-positive definiteness.

mchol(G)
Cholesky decomposition, especially to have PD matrix.

ltqnorm(p)
Lower tail quantile for standard normal distribution. it returns the Z satisfying Pr{X < Z} = P

ranking(x)  
Ranking of a vector


Pattern Generation help

drawPattern(IM)
Draws one pattern according to colormap mycamp, and plots it beautifully, X and Pat are assumed global variables.

[distance,index] = find_closest_pattern_index(randPattern_inR,n)
Finds closest n patterns close to randPattern_inR and then gives their indexes with the corresponding distances.

find_corresponding_points(ptsA,ptsB)
Finds corresponding control points of B in A and outputs them in a matrix  with the dimension equal to the number of points in B.

find_skel_ends
Finds skeleton end points. For example of usage:
I = imread(filename);
input_skeleton_image = bwmorph(I,'thin',inf);
terminating_pts = find_skel_ends(input_skeleton_image);
figure,imshow(input_skeleton_image);
hold on; 	plot(terminating_pts(:,1),terminating_pts(:,2),'r*');
find_skel_intersection
Finds skeleton intersection points

morphimage(orig,morphfield)
Morphs the image with a morphfield, used in pattern_morph function. 

I = skeleton_CP_warp(i,j,w)
Skeleton warping on pattern I and j , using a factor of w. uses warpImage function too. As an example:
[distance,index] = find_closest_pattern_index(randPattern_inR,2);
I = skeleton_CP_warp( index(1), index(2),  1/(distance(1)*sum(1./distance)) );
Pattern = I>0.5; % issue with matlab (?)
figure;imshow(Pattern);

pattern_morph(i,morphIntensity)
Morph a pattern. An example is given here:
[distance,index] = find_closest_pattern_index(randPattern_inR,1);
Pattern = pattern_morph(index(1),1.5);
figure;imshow(Pattern);

[mappedX, mapping] = pcanalysis(X, no_dims)
Performs pca analyss. Used in KL decomposition as given in example below (Not Used anymore):
[distance,index] = find_closest_pattern_index(randPattern_inR,1);
[mappedX, mapping] = pcanalysis(X, Pat^2);
reconst=mappedX*mapping.M';

[I,threshold] = Radon_warp(distance,index)
Radon warping indexed X’s according to distance and giving the threshold.For example:
[distance,index] = find_closest_pattern_index(randPattern_inR,3);
[I,threshold] = Radon_warp(distance,index);
Pattern = I>threshold;
figure;imshow(Pattern);

warpImage(Image, originalMarks, desiredMarks)
Warps Image so that orinal marked points are at the desired marked points and all other points 
are interpolated in correctly.

warpitfun(pink,get,Zp,Zg)
The function a very simple Thin-Plane spline warping of the two images. (Not used as it does not do interpolation)


Other unknown functions

dijk(A,s,t)

dijkstra( G , S )

sortclasses(given_inputs, separation_distance, connectivity)

statgetargs
Some function used somewhere to check vararging parameters for statistics function.

stprpath(toolboxroot)
adds path of toolbox to the Matlab start up file.










