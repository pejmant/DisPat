avgReproduction = 0;

dm = floor(sqrt(clus));
if dm ~= sqrt(clus)
    dm = dm+1;
end   

temp = zeros(disDim2,Pat^2);
prototype = zeros(clus,Pat^2);
for i=1:disDim2
    DT=bwdist(reshape(X(i,:),Pat,Pat));
    if DT==Inf
        DT = ones(Pat,Pat);
    end
    temp(i,:)=reshape(1-DT./max(max(DT)),1,Pat^2);
end

figure;
for i = 1:clus
    clust = find(idx==i);
    n=size(clust,1);
    subplot(dm,dm,i);
    prototype(i,:)=sum(temp(clust,:),1)./n;
    imshow(reshape(prototype(i,:),Pat,Pat),[0.3 1],'InitialMagnification','fit');

    reProd = sum(abs(2.*prototype(i,:)-1))/Pat^2;
    avgReproduction = avgReproduction + reProd*n/disDim2;
    
    xlabel(num2str(reProd));
    title(['\bfCluster ',num2str(i)]);
end

fprintf('\nweighted sum of sharpnesses  of  patterns (proximity)= %2.6f\n',avgReproduction);

