function [distance,index] = find_closest_pattern_index(randPattern_inR,n)
    global Y;
    d = zeros(size(Y,1),1);
    for j=1:size(Y,1)
        d(j,1) = pdist(vertcat(Y(j,:),randPattern_inR));
    end
   [distance,index] = sort(d,'ascend');
   index = index(1:n,1);
   distance = distance(1:n,1);
end