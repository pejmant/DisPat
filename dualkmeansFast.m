function f = dualkmeansFast(K,N,idx)

%function [f,d] = dualkmeans(K,N)
%
% Performs dual K-means for Nr samples specified by the kernel K
%
%INPUTS
% K = the kernel matrix
% N = the number of clusters desired
% sigma =
% idx =
%
%OUTPUTS
% f = the cluster allocation vector
% d = the distances of the samples to their respective cluster centroids
%
%
%For more info, see www.kernel-methods.net


% original kernel matrix stored in variable K
% clustering given by a Nr x N binary matrix A
% and cluster allocation function f
% d gives the distances to cluster centroids

Nr=size(K,1);
A = zeros(Nr,N);
%f = ceil(rand(Nr,1)* N); 

% If no initial points are specified - random
if nargin == 2,
    f = ceil(rand(Nr,1)* N);
    for i=1:Nr
      A(i,f(i)) = 1;
    end
end

%Use initial points to construct A
if nargin == 3, 
	if size(idx,2) == N   %idx contains the indexes of the centroids among the Nr Points
		dist = zeros(1,N);
		for i =1:Nr
			for j = 1:N
				dist(j)=featuredist(i,idx(j),K); 
			end
			[B,IX]=min(dist);
			A(i,IX)=1;  % A(i,j) = 1 if the realization i belongs to cluster j, 0 otherwise
			f(i) = IX;  % f is a vector of the indexes of clusters for each points
		end
	else
		if size(idx,1) == Nr
			for i =1:Nr
				A(i,idx(i)) = 1;
			end
			f=idx;
		else
			%'Error - Number of Initial Centroids different from number of clusters!';
		end
	end
end


change = 1;
iter = 0;
while change == 1 && iter < 100 
  change = 0;
  E = A * diag(1./sum(A));
  KE = K*E;
  Z = bsxfun(@minus,ones(Nr,1)* diag(E'*KE)', 2*KE);
  [d, ff] = min(Z, [], 2);
  for i=1:Nr
   if f(i) ~= ff(i)
    A(i,ff(i)) = 1;
    A(i, f(i)) = 0;
    change = 1;
   end
  end
  f = ff;
  iter = iter + 1;
end

