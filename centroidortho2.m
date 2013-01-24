function idxclusters = centroidortho2(Y,k)

% This function 

% Input Parameter:
% Y = Matrix of k orthonormal eigenvectors (by colomn)
% k = Number of cluster - number of "orthogonal" points to be constructed 

% Output Paramater
% idxclusters: vector containing the index of points that are the most
% orthogonal

[mi m] = min(Y(:,1));
[Mi M] = max(Y(:,1));

% Start by choosing the extremes according to the first axe
idxclusters = [m M]';

for i = 3:k
    Ybis=Y;
    Ybis(idxclusters,:)=[];
    [ps t] = min(sum(Y(idxclusters,:)*Ybis'));  % Minimization of the dot product between points
    idxclusters = vertcat(idxclusters,t);
end