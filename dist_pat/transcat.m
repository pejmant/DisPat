function realization_C = transcat(realization, out, hardData, par, omega)

% This function will apply TRANSCAT to categorical training images in order
% to apply correction on proportions.       -- Mehrdad Honarkhah --

lim  = (par.Pat  - 1)/2;
limz = (par.Patz - 1)/2;
realization_C = realization;


% fixing the limits for categorical variables (more than 2)
if sum(out(:)) ~= sum(round(out(:)))
    out         = round(out.*(length(unique(out(:))) -1));
    realization = round(realization.*(length(unique(out(:))) -1));
end



% -----------------------------------------------------------------
% find the proportions initially in realization and training image
% -----------------------------------------------------------------
% calculating proportions seen in TI for all categorical variables
proportions.TI = zeros(1,max(out(:))+1);
for i = 1:max(out(:))+1
    proportions.TI(i) = sum(out(:) == i-1)/(par.Dim^2*par.Dimz);
end

par.Dim  = size(realization,1);
par.Dimz = size(realization,3);

% calculating proportions seen in realization for all categorical variables 
proportions.Re = zeros(size(proportions.TI));
for i = 1:length(proportions.TI)
    proportions.Re(i) = sum(realization(:) == i-1)/(par.Dim^2*par.Dimz);
end





% -----------------------------------------------------------------
% design the filter
% -----------------------------------------------------------------
filter = ones(par.Pat,par.Pat,par.Patz);

mainAxis  = [1:2:par.Pat , par.Pat-2:-2:1];
if par.Patz == 1
    mainAxisz = par.Pat;
else
    if par.Patz < par.Pat
        taxis = round(linspace(1,par.Pat,(par.Patz+1)/2));
        mainAxisz = [taxis,taxis(end-1:-1:1)];
    else
        mainAxisz = [1:2:par.Patz, par.Patz-2:-2:1];
    end
end

filter(1:par.Pat , lim+1     , limz+1    ) = mainAxis' ;
filter(lim+1     , 1:par.Pat , limz+1    ) = mainAxis  ;
filter(lim+1     , lim+1     , 1:par.Patz) = mainAxisz ;
Vc = max(par.Pat,par.Patz);
for i =1:par.Pat
    for j=1:par.Pat
        for k =1:par.Patz
            if i==lim +1 || i==1 || i==par.Pat,  continue; end;
            if j==lim +1 || j==1 || j==par.Pat,  continue; end;
            if par.Patz > 1
                if k==limz+1 || k==1 || k==par.Patz, continue; end;
            end
            filter(i,j,k) = round(Vc/(1+sqrt((i-lim-1)^2+(j-lim-1)^2+(k-limz-1)^2)));
        end
    end
end




% -----------------------------------------------------------------
% make hard data conditioning factor
% -----------------------------------------------------------------
hardData(~isnan(hardData)) = 10;
hardData(isnan(hardData))  = 1;


% -----------------------------------------------------------------
% change of proportions to be applied to realization
% -----------------------------------------------------------------

% omega = (new concept not to allow big proportion change, i.e. max change = 10%)
propChange = abs(proportions.TI - proportions.Re)./proportions.TI;
for i = 1:length(proportions.TI)
    if propChange(i) > omega
        propChange(i) = 1 + 0.1*sign(proportions.Re(i)-proportions.TI(i));
    else
        propChange(i) = proportions.TI(i)./proportions.Re(i);
    end
end


% to temporarily hold the proportions in each window
tempProp = zeros(size(proportions.TI));



% -----------------------------------------------------------------
% random walk through the nodes
% -----------------------------------------------------------------
randm = randperm(par.Dim^2*par.Dimz);
% randm = 1:par.Dim^2*par.Dimz;

for r = 1:par.Dim^2*par.Dimz
   [i,j,k] = ind2sub([par.Dim, par.Dim, par.Dimz], randm(r));

   % calculation to consider for boundary nodes
   rx1 = i-lim; ry1 = j-lim; rz1 = k-limz;
   rx2 = i+lim; ry2 = j+lim; rz2 = k+limz;
   fx1=1; fy1=1; fz1=1;
   fx2=par.Pat; fy2=par.Pat; fz2=par.Patz;
   
   % left boundaries
   if i-lim < 1
       rx1 = 1;
       rx2 = par.Pat - (lim - i) -1;
       fx1 = lim-i+2;
   end
   if j-lim < 1
       ry1 = 1;
       ry2 = par.Pat - (lim - j) -1;
       fy1 = lim-j+2;
   end
   if k-limz < 1
       rz1 = 1;
       rz2 = par.Patz- (limz- k) -1;
       fz1 = limz-k+2;
   end       
   
   % right boundaries
   if i+lim > par.Dim
       rx1 = i+1+lim-par.Pat;
       rx2 = par.Dim;
       fx2 = par.Pat - (i+lim -par.Dim );
   end
   if j+lim > par.Dim
       ry1 = j+1+lim-par.Pat;
       ry2 = par.Dim;
       fy2 = par.Pat - (j+lim -par.Dim );
   end
   if k+limz > par.Dimz
       rz1 = k+1+limz-par.Patz;
       rz2 = par.Dimz;
       fz2 = par.Patz- (k+limz-par.Dimz);
   end 

   % make necessary changes for the limits (accorording to boundary or not)
   temp       = realization(rx1:rx2 , ry1:ry2 , rz1:rz2);
   filterTemp = filter     (fx1:fx2 , fy1:fy2 , fz1:fz2);
   hardTemp   = hardData   (rx1:rx2 , ry1:ry2 , rz1:rz2);
   
   % apply filter
   for l = 1:length(proportions.TI)
       B = (temp==l-1);
       B = +B;
       tempProp(l) = mean(B(:).*filterTemp(:).*hardTemp(:));
   end
   tempProp = tempProp./sum(tempProp(:));
   tempProp = tempProp.*propChange;
   [dummy, indexProp] = max(tempProp);
   realization_C(i,j,k) = indexProp - 1;
   

end


% fixing back the max and min to be in:[0,1] range
realization_C = realization_C./(length(unique(out(:))) -1);

end