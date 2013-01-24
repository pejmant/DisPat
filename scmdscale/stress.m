function y = stress(D1,D2)
%y = stress(D1,D2);
%D1 is the distance matrix; for example: D1 = pdist(X);
%D2 is the compared distance matrix or disparities
%D1 and D2 is computed by pdist function

if numel(D1)~=numel(D2)
	error('size of D1 and D2 must be the same');
end

maxD1 = max(max(D1));
maxD2 = max(max(D2));

s = maxD2/maxD1;
D2 = D2/s;

y = sqrt(sum((D1-D2).^2)/sum(D1.^2));