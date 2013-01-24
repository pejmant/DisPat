disDim2=size(X,1);
disMat = zeros(disDim2);
disp('Which dissimilarity distance method ? ') ;
disp(' 1 . Euclidean distance ');
disp(' 2 . Hausdorff distance ');
disp(' 3 . Modified Hausdorff distance ');
disp(' 4 . Hough Transform dissimilarity distance');
disp(' 5 . Cosine distance');
disp(' 6 . Radon Transform dissimilarity distance');
disp(' 7 . distance Transform ');
disp(' 8 . distance Transform (proximity)');
way = input(' ');


switch way
    case 8
        temp=zeros(disDim2,Pat^2);
        for i=1:disDim2
            DT=bwdist(reshape(X(i,:),Pat,Pat));
            if DT==Inf
                DT = ones(1,Pat^2);
            end
            temp(i,:)=reshape(1-DT./max(max(DT)),1,Pat^2);
        end
        disMat=squareform(pdist(temp));
        
    case 7
        temp=zeros(disDim2,Pat^2);
        for i=1:disDim2
            temp(i,:)=reshape(bwdist(reshape(X(i,:),Pat,Pat)),1,Pat^2);
            if temp(i,:)==Inf
                temp(i,:)= sqrt(2)*(Pat-1).*ones(1,Pat^2);
            end
        end
        disMat=squareform(pdist(temp));
    
    case 6
        fprintf('\n Percentage Completed:    ');
        theta=0:359;
        for i=1:disDim2
            perc=round(100*i/disDim2);
            fprintf('\b\b\b\b%3d%%', perc);
            temp(i,:)=reshape(radon(reshape(X(i,:),Pat,Pat),theta),1,[]);
        end
        fprintf('\n>Calculating distances...');
        disMat=squareform(pdist(temp));
        
    case 5
        fprintf('>Calculating distances...');
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            for j=i+1:disDim2
                disMat(i,j) = calcDissimilarityCosine(i,j,disDim2,Pat);
                disMat(j,i) = disMat(i,j);
            end
        end
        
     case 4
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            perc=round(100*i/disDim2);
            fprintf('\b\b\b\b%3d%%', perc);
            temp(i,:)=reshape(hough(reshape(X(i,:),Pat,Pat)),1,[]);
        end
        fprintf('\n>Calculating distances...');
        disMat=squareform(pdist(temp));

     case 3
        fprintf('>Calculating distances...');
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            for j=i+1:disDim2
                disMat(i,j) = calcDissimilarityModified(i,j,disDim2,Pat);
                disMat(j,i) = disMat(i,j);
            end
        end    

     case 2
        fprintf('>Calculating distances...');
        fprintf('\n Percentage Completed:    ');
        for i=1:disDim2
            for j=i+1:disDim2
                disMat(i,j) = calcDissimilarity(i,j,disDim2,Pat);
                disMat(j,i) = disMat(i,j);
            end
        end
        
        
    case 1       
         disMat=squareform(pdist(X));
%         disMat=squareform(pdist(X),'seuclidean');
%         disMat=squareform(pdist(X),'mahalanobis');
%         disMat=squareform(pdist(X),'cityblock');
%         disMat=squareform(pdist(X),'minkowski',3);
%         disMat=squareform(pdist(X),'cosine');
%         disMat=squareform(pdist(X),'correlation');
%         disMat=squareform(pdist(X),'spearman');
%         disMat=squareform(pdist(X),'hamming');
%         disMat=squareform(pdist(X),'jaccard');
%         disMat=squareform(pdist(X),'chebychev');
        
        
end