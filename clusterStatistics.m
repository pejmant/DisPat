
function [cluster]=clusterStatistics(Z,i,idx)
% statistics of the clusters in kernel space
% mean and std and corr are calculated using
% the kernel model. Kernel PCA is used to project
% data into a Pat^2 dimensions (No Error).
    
    clust_show = find(idx==i);
    cluster.mean=mean(Z(clust_show,:));
    cluster.corr =corrcoef(Z(clust_show,:));
    cluster.std =std(Z(clust_show,:));
    
end
