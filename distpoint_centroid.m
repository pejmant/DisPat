function Dist = distpoint_centroid(ff,K,N)

% Function calculating the distance between the centroid and the points of the cluster

ell = size(K,1);
Dist = zeros(ell,N);

for i = 1:ell
    for j = 1:N
    	dim=size(find(ff==j),1);
    	if dim == 0
    		dim = 1;
    	end
        Dist(i,j) = K(i,i) - 2/dim*sum(K(i,find(ff==j))) + sum(sum(K(find(ff==j),find(ff==j))))/(dim^2);
    end
end