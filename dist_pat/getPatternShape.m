function [wx, wy, wz] = getPatternShape(node, par)

% This function gets the wx and wy  and wz which reflects the array where the
% pattern should be read from or written to. In other words, this array of
% wx and wy points to the location where the multiple-grid pattern can be
% created from.

Pat           = par.Pat;
Patz          = par.Patz;
m1            = par. m1;
multipleGrid  = par.multipleGrid;
node = node';

% to account for created borders around realization
borders = (Pat  - 1)*2^(multipleGrid-1)/2;   
bordersz = (Patz - 1)*2^(multipleGrid-1)/2; 

% to make the grid be at center location of the template
r  = (Pat  + 1)/2;    
rz = (Patz + 1)/2;

i  = m1.*(0:Pat -1);
iz = m1.*(0:Patz-1);


wx = bsxfun(@times, ones(size(node,1), Pat),  1 + m1*(node(:,1) -r) + borders);
wx = bsxfun(@plus, wx, i);
 
wy = bsxfun(@times, ones(size(node,1), Pat),  1 + m1*(node(:,2) -r) + borders);
wy = bsxfun(@plus, wy, i);

wz = bsxfun(@times, ones(size(node,1), Patz), 1 + m1*(node(:,3)-rz) + bordersz);
wz = bsxfun(@plus, wz, iz);


end