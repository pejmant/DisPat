% Compute the distance between two points in the feature space
%||Phi(x) - Phi(z)||2 = k(x,x) - 2k(x,z) + k(z,z)



function dist = featuredist(i,j,K)
	dist = K(i,i) + K(j,j) -2*K(i,j);
	
	