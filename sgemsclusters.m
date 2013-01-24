function sgemsclusters
    Dim = 51;
    Pat = 9;
    [datain, colnames, line1]=loadgeoeas('prot51.dat');
    proto = geoeas2matlab(datain,[Dim Dim]);
    cluss = max(max(proto))+1;
    
    for i=1:Dim*Dim
        if proto(i)==-9999
        proto(i)=-1;
        end
    end
    imagesc(proto);
    
    
    image=imread('ti.bmp');
    image=imresize(image,[Dim Dim]);
    image=im2bw(image,0.75);
    out=1-image;

    disDim = Dim - Pat + 1;
    X=[];
    k=1;
    for i=1:disDim
        for j=1:disDim
            X(k,:)=reshape(out(i:i+Pat-1,j:j+Pat-1),1,Pat^2);
            k=k+1;
        end
    end
    
    RegionProto = proto(floor(Pat/2)+1:Dim-floor(Pat/2),floor(Pat/2)+1:Dim-floor(Pat/2));
    figure; imagesc(RegionProto);
    idx = reshape(RegionProto',disDim*disDim,1);
    avgReproduct = 0;
    avgReproduct1 = 0;
    figure;
    dm = floor(sqrt(cluss));
    if dm ~= sqrt(cluss)
        dm = dm+1;
    end        
    for i = 1:cluss
        clust = find(idx==i-1);
        n=size(clust,1);
        subplot(dm,dm,i);
        imshow(reshape(sum(X(clust,:),1),Pat,Pat),[0.3*n n]);
        rePro = sum(abs(2.*(sum(X(clust,:),1)./n)-1))/Pat^2;
        avgReproduct = avgReproduct + rePro;
        avgReproduct1 = avgReproduct1 + rePro*n/1849;
        xlabel(num2str(rePro));
        title(['\bfCluster ',num2str(i)]);
    end
    fprintf('\nAverage exactness of patterns in clusters = %2.6f\n',avgReproduct/cluss);
    fprintf('\nweighted sum of sharpnesses  of  patterns= %2.6f\n',avgReproduct1);
    
    
end
