function [I,threshold] = Radon_warp1(distance,index,radonX)
    global X;global Pat;
    theta=0:359;
    distance=(1./distance).^2;
    AR = distance(1).*squeeze(radonX(index(1),:,:));
    meanI = distance(1)*mean(X(index(1),:));
    for i=2:length(distance)
        AR    = AR    + distance(i).*squeeze(radonX(index(i),:,:));
        meanI = meanI + distance(i)*mean(X(index(i),:));
    end
    AR = AR./sum(distance);
    meanI = meanI./sum(distance);
    I = iradon(AR, theta, 'linear','Ram-Lak', 1, Pat);
    threshold = quantile(reshape(I,1,[]),1-meanI);
end