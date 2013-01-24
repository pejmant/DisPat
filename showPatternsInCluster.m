function showPatternsInCluster(i,idx)
    global X;global Pat;
    clust_show = find(idx==i);
    n=size(clust_show,1);
    dmn = ceil(sqrt(n));
    for j=1:n
       subplot(dmn,dmn,j);
       imshow(reshape(X(clust_show(j),:),Pat,Pat));
       title(['Cluster ',num2str(clust_show(j))]);
    end
end
