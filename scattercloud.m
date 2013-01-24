function h = scattercloud(Y)
%SCATTERCLOUD display density of scatter data
%   SCATTERCLOUD(X,Y) creates a scatterplot of X and Y, displayed over a
%   surface representing the smoothed density of the points.  The density is
%   determined with a 2D histogram, using 25 equally spaced bins in both
%   directions.
%   SCATTERCLOUD(X,Y,N) uses N equally spaced bins.
%   SCATTERCLOUD(X,Y,N,L) uses L as a parameter to the smoothing algorithm.
%    Defaults to 1.  Larger values of L lead to a smoother density, but a
%    worse fit to the original data.
%   SCATTERCLOUD(X,Y,N,L,CLM) uses CLM as the color/linestyle/marker for
%    the scatter plot.  Defaults to 'k+'.
%   SCATTERCLOUD(X,Y,N,L,CLM,CMAP) uses CMAP as the figure's colormap.  The
%    default is 'flipud(gray(256))'.
%   H = SCATTERCLOUD(...) returns the handles for the surface and line
%    objects created.
%
%   Example:
%
%     scattercloud(1:100 + randn(1,100), sin(1:100) + randn(1,100),...
%                  50,.5,'rx',jet(256))

x=Y(:,1);y=Y(:,2);z=Y(:,3);
x = x(:);
y = y(:);
z=z(:);
cmap = flipud(gray(256));
clm = 'k+';
l = 1;  
n = 20;

% min/max of x and y
minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);
minZ = min(z);
maxZ = max(z);

% edge locations
xEdges = linspace(minX,maxX,n);
yEdges = linspace(minY,maxY,n);
zEdges = linspace(minZ,maxZ,n);

% shift edges
xDiff = xEdges(2) - xEdges(1);
yDiff = yEdges(2) - yEdges(1);
zDiff = zEdges(2) - zEdges(1);
xEdges = [-Inf, xEdges(2:end) - xDiff/2, Inf];
yEdges = [-Inf, yEdges(2:end) - yDiff/2, Inf];
zEdges = [-Inf, zEdges(2:end) - zDiff/2, Inf];
% number of edges
numX = numel(xEdges);
numY = numel(yEdges);
numZ = numel(zEdges);
% hold counts
C = zeros(numY,numX,numZ);

% do counts
for i = 1:numY-1
    for j = 1:numX-1
        for k = 1:numZ-1
            C(i,j,k) = length(find(x >= xEdges(j) & x < xEdges(j+1) &...
                                   y >= yEdges(i) & y < yEdges(i+1) &...
                                   z >= zEdges(k) & z < zEdges(k+1)));
        end
    end
end

% get rid of Infs from the edges
xEdges = [xEdges(2) - xDiff,xEdges(2:end-1), xEdges(end-1) + xDiff];
yEdges = [yEdges(2) - yDiff,yEdges(2:end-1), yEdges(end-1) + yDiff];
zEdges = [zEdges(2) - zDiff,zEdges(2:end-1), zEdges(end-1) + zDiff];
% smooth the density data, in both directions.
  C = (0+2*smooth3(C,'gaussian',[3,3,3]))/2 ;

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'ZMinorTick','on',...
    'ZColor',[0.502 0.502 0.502],...
    'YMinorTick','on',...
    'YColor',[0.502 0.502 0.502],...
    'XMinorTick','on',...
    'XColor',[0.502 0.502 0.502],...
    'Projection','perspective',...
    'FontWeight','demi',...
    'FontSize',8,...
    'FontName','Garamond');

view([20 42]);
box('on');
grid('on');
hold('all');

p = plot3(x,y,z,'k.');
axis('tight');
hold on;
h = slice(xEdges,yEdges,zEdges,C,xEdges,yEdges,zEdges);
alpha('color')
set(h,'EdgeColor','none','FaceColor','interp',...
	'FaceAlpha','interp')
% colormap(cmap);

grid('on');
box('on');





