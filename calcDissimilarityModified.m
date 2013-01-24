function disMat = calcDissimilarityModified(i,j,disDim2,Pat)
global X;

%% find the two corresponding patterns A and B
A = reshape(X(i,:),Pat,Pat);
B = reshape(X(j,:),Pat,Pat);

% Progress Display
perc=round(100*((i-1)*disDim2-i*(i-1)/2+j-i)/(disDim2*(disDim2-1)/2));
fprintf('\b\b\b\b%3d%%', perc);

%% Calculate dissimilarity by Modified Hausdorff distance

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