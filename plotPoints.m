
figure;
ptsymb = {'bs','r^','md','go','c+','y*','kv','b>','r<','mp','gh','c.'};
ptsymb=repmat(ptsymb,1,8); PC=0;



for i = 1:clus
    clust = find(idx==i);
    plot3(Z(clust,1),Z(clust,2),Z(clust,3),ptsymb{i});
    hold on
end

title('KPCA results in 3-D kernel space with centers');
grid on;
