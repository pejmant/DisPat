function disMat = haussdorff(A,B)


% Calculate dissimilarity by Modified Hausdorff distance

    [r,c] = find(A);
    pointA = [r c];
    [r,c] = find(B);
    pointB = [r c];
    disMat = (max( compute_dist(pointA,pointB), compute_dist(pointB,pointA) ));


            
function[ dist ] = compute_dist( A, B )

m = size(A,1);
n = size(B,1);
dist = zeros(m,1);
for k = 1 : m
    C = ones(n,1) * A(k,:);
    D = (C-B).^2;
    dist(k) = min(sum(D,2));
end
dist = sum(sqrt(dist))/m;