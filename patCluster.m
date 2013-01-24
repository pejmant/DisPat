% Choosing between clustering methods
disp('Clustering Method ? ') ;
disp('Kernel K-means analysis ');
way2 = 3;
MDS1=MDS;

clus = input('How many clusters ? ') ;

sigma =0.2*max(max(pdist(Y)));

% Spectral Clustering to find initial clusters
[idx,XX,totsumD] = SpectralClustNg(Y,sigma,clus);

% Apply Radial Basis Function Kernel
K = kernel(Y','rbf',sigma);

% Define the dimension of the projections on the Feature space
% e = eig(K);
%new_dim=length(find(e>1));
new_dim=3;
model= kpca(Y', struct('ker','rbf','arg',sigma,'new_dim',new_dim));
Z = kernelproj(Y',model)';


% Kernel K-Means algorithm with initial points from Spectral Clustering.
[idx,f,A] = dualkmeans3(K,clus,sigma,idx);


MDS=-1;