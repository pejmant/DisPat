function result = kmedoid(X,c)

N = size(X,1);

%initialization
index=randperm(N);
v=X(index(1:c),:);v = v + 1e-10;
v0=X(index(1:c)+1,:);v0 = v0 - 1e-10;
% dist=zeros(size(X,1),c);

vIndex = zeros(c,1);

while prod(max(abs(v - v0)))
    
        v0 = v;
        
        %Calculating the distances
        dist = slmetric_pw(X', v', 'sqdist');
        
        %Assigning clusters
        [m,label] = min(dist,[],2);
        label = label';

      
        %Calculating cluster centers
        for i = 1:c
            index=find(label == i);
            if ~isempty(index)  
                vtemp       = mean(X(index,:));
                [dummy, inx]= min(sum(bsxfun(@minus, X, vtemp).^2,2));
                v(i,:)      = X(inx,:);
                vIndex(i,1) = inx;
            else 
                v(i,:)=X(round(rand*N-1),:);
            end   
        end
end

%results
result.v     = v;
result.index = vIndex;

end
