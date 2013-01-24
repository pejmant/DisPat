if way2==2
    dm = floor(sqrt(numClust));
    if dm ~= sqrt(numClust)
        dm = dm+1;
    end        
    for i = 1:numClust
        clust = clustMembsCell{i};
        n=size(clust,2);
        subplot(dm,dm,i);
        imshow(reshape(sum(X(clust,:),1),Pat,Pat),[0.3*n n]);
        title(['Cluster ',num2str(i)]);
    end
    MDS=MDS1;
end
% If K means
if way2==1
    calcClusters
end
% Kernel K-means
if way2==3
    calcClusters
    
%     figure;
%     for i = 1:clus
%         subplot(dm,dm,i);
%         imshow(reshape(X(clusterCenterIndex(i),:),Pat,Pat));
%         title(['Cluster ',num2str(i)]);
%     end

    MDS=MDS1;
end