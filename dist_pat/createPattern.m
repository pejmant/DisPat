function [Pattern, cluster] = createPattern(n, clusterModel, Z, X, Y, idx, Pat, cluster,radonX)

global Y;
global X;
global Pat;
%% Find Random Point in a Cluster

% Finding a Point in cluster using the model
cluster = clusterStatistics(Z, n, idx);

randPattern = lhs_iman_n(cluster.mean, cluster.std, cluster.corr);

randPattern_inR = Yrec(randPattern', clusterModel)';



%% Pattern warping
% using Radon Algorithm
[distance,index] = find_closest_pattern_index(randPattern_inR,3);
[I,threshold] = Radon_warp1(distance,index,radonX);
Pattern = I>threshold;
% figure;drawPattern(Pattern);