function radonX = calculateRadonX(X, Pat)

theta=0:359;

fprintf('Calculating Radon of Patterns .................');

[sz1,sz2] = size(radon(reshape(X(1,:),Pat,Pat),theta));
radonX = zeros(size(X,1),sz1,sz2);

for i=1:size(X,1)
    radonX(i,:,:) = radon(reshape(X(i,:),Pat,Pat),theta);
end
  
fprintf(' Done!\n');

end