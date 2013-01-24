figure;
Prototype=zeros(Dim,Dim);
for i = 1:clus
    clust = find(idx==i);
    Iclust = ceil(clust./disDim);
    Jclust = mod(clust,disDim)+1;
    for j = 1:length(Iclust)
        Prototype(1+m.*(Iclust(j)-1)+(Pat-1)*m1/2,1+m.*(Jclust(j)-1)+(Pat-1)*m1/2)=i;
    end
end
imagesc(Prototype);
axis square
