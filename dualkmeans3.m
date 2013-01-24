function [f,d,A,b] = dualkmeans3(K,N,sigma,idx)

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
if nargin == 3,
	f = zeros(Nr,1);
	while(size(unique(f),1) < N)
		f = ceil(rand(Nr,1)* N); % f(i) is a vector containing the index of where the 1 is in each row of A.
	end
	for i=1:Nr
	  A(i,f(i)) = 1;
    end
end

%Use initial points to construct A
if nargin == 4, 
		
	if size(idx,2) == N   %idx contains the indexes of the centroids among the Nr realizations
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
			for i =1:size(K,1)
				A(i,idx(i)) = 1;
			end
			f=idx;
		else
			'Error - Number of Initial Centroids different from number of clusters!'
		end
	end
end


[a b] = max(A,[],2);

change = 1;

while change == 1

  changevect=[];
  l=1;
  emptycluster = [];

  Z = distpoint_centroid(f,K,N);
  [d, ff] = min(Z, [], 2);

  if size(unique(ff)) < N % emptycluster = which point not represented
      for i = 1:N
            if any(i == ff) == 0
                    emptycluster(l) = i;
                    l=l+1;
            end
      end
      [dd, fff] =  min(Z, [], 1); 
      ff(fff(emptycluster)) = emptycluster;
  end

  for i=1:Nr
   if (f(i) ~= ff(i)) 
       A(i,ff(i)) = 1;
       A(i, f(i)) = 0;
       changevect(i)=1;
   else
       changevect(i)=0;
   end
  end

  change = any(changevect);
  f = ff;
end
