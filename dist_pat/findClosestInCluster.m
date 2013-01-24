function [Pattern, patternIdx] = findClosestInCluster(dataEvent, X, index, par, radonX, wieghtEvent)

% This function finds the closest pattern to the dataevent amongst all the
% pattern in a specific cluster


Xtemp = X(index,:);


difference = sum(bsxfun(@times,abs(bsxfun(@minus, Xtemp, dataEvent)),wieghtEvent), 2);
[dummy, idxNumber] = min(difference);


if par.bDataEventOptimization && dummy ~= 0
% not for 3D case
        % Finding the second minimum one
        % ------------------------------
        [B,IY] = sort(difference);
        difference(IY(B == B(1))) = inf;
        [dummy2, idxNumber2] = min(difference);

        % make the minimization variables and minimize
        Pat = par.Pat;
        reshapedDataEvent = reshape(dataEvent, Pat, Pat);
        X1 = X(idxNumber,:);
        X2 = X(idxNumber2,:);
        radonX1 = squeeze(radonX(idxNumber,:,:));
        radonX2 = squeeze(radonX(idxNumber2,:,:));
        alpha0 = 0.5;
        f = @(alpha)objectiveFunc(alpha, reshapedDataEvent, X1, X2, radonX1, radonX2, Pat);

        % finding the minimum 
        options = optimset('Display','off');
        [minAlpha, fval]= fminbnd(f, 0, 1, options);

        % check if it is better than original pattern
        if fval < dummy
                % Assigning the Pattern to be pasted
                theta=0:359;
                AR = minAlpha*radonX1 + (1-minAlpha)*radonX2;
                meanI = minAlpha*mean(X1) + (1-minAlpha)*mean(X2);
                I = iradon(AR, theta, 'linear','Ram-Lak', 1, par.Pat);
                threshold = quantile(reshape(I,1,[]),1-meanI);
                Pattern = I>threshold;
                Pattern = reshape(Pattern, 1, par.Pat^2);
        else
                Pattern = Xtemp(idxNumber,:);
        end

        % find the range of the dual template
        patternIdx = findRangeDualTemplate(index(idxNumber), par);

else
        % Assigning the Pattern to be pasted
        Pattern = Xtemp(idxNumber,:);

        % find the range of the dual template
        patternIdx = [];
        if par.bUseDualTemplate
            patternIdx = findRangeDualTemplate(index(idxNumber), par);
        end
end



end






function output = objectiveFunc(alpha, reshapedDataEvent, X1, X2, radonX1, radonX2, Pat)

    theta=0:359;
    AR = alpha*radonX1 + (1-alpha)*radonX2;
    meanI = alpha*mean(X1) + (1-alpha)*mean(X2);
    I = iradon(AR, theta, 'linear','Ram-Lak', 1, Pat);
    threshold = quantile(reshape(I,1,[]),1-meanI);
    Pattern = I>threshold;


    output = sum(abs(Pattern - reshapedDataEvent), 2);

end


