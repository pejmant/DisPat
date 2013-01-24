function T = SpectralClustNg(K,N)

% Performs spectral K-means for ell samples specified by the kernel K
%
%INPUTS
% K = the kernel matrix
% N = the number of clusters desired
%
%OUTPUTS
% T = vector (Nr*1) containing, for each point, the index of the cluster he belongs to  
% X = Matrix (Nr*N) containing the eigenvectors of L 
% totsumD = Total sum of values of distances from each point to every centroid.

% Reference: 
% - Dhillon, I., Guan Y and Kulis B., Kernel k-means, Spectral Clustering and Normalized Cuts, 
% KKD'04, August 22-25, 2004, Seattle, Washinton, USA



invSqrtD = diag(1./sqrt(sum(K,1)));
L = invSqrtD*K*invSqrtD;

% Find the N largest eigenvectors of L
% [V,D] = eig(L);
% X = V(:,1:N);
% Method Mehrdad for faster eig
[X,D] = symeig(L, N);


% Normalization of X
Z = sqrt(sum(X.^2, 2));
% X = X ./ Z(:, ones(size(X,2), 1));
X = bsxfun(@rdivide,X,Z);

% Try to find N points that are orthogonal to initialize KKM and avoid
% randomness of the results
idxclusters = centroidortho2(X,N);

% Perform traditionnal k-means algorithm using previous initial points
[T a sumD D totsumD] = kmeans2(X,N,'Start',X(idxclusters,:),'Maxiter',500,'EmptyAction','singleton','Display','final');

