function disMat = calcDissimilarityCosine(i,j,disDim2,Pat)
global X;

% Progress Display
perc=round(100*((i-1)*disDim2-i*(i-1)/2+j-i)/(disDim2*(disDim2-1)/2));
fprintf('\b\b\b\b%3d%%', perc);

%% Calculate dissimilarity by a measure of Cosine distance
disMat = acos(sum(X(i,:).*X(j,:))/(sqrt(sum(X(i,:)))*sqrt(sum(X(j,:)))));
