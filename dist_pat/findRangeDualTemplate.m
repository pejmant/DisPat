function patternIdx = findRangeDualTemplate(index, par)


m= par.m;
m1 = par.m1;
Dim = par.Dim;
Dimz= par.Dimz;
Pat = par.Pat;
Patz= par.Patz;

disDim    = ceil((Dim - (1+(Pat -1)*m1) + 1)/m);
disDimz   = ceil((Dimz -(1+(Patz-1)*m1) + 1)/m);
[i, j, k] = ind2sub([disDim, disDim, disDimz], index);
patternIdx.x = 1+m*(i-1):1+m*(i-1)+(Pat -1)*m1;
patternIdx.y = 1+m*(j-1):1+m*(j-1)+(Pat -1)*m1;
patternIdx.z = 1+m*(k-1):1+m*(k-1)+(Patz-1)*m1;

end