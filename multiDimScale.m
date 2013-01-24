fprintf('\n');
disp('MultiDimensional Scaling ') ;
disp('Classical  MDS');

[Y1,e] = cmdscale(double(disMat));


% Plot EigenValues
    figure;
    subplot(2,1,1);
    plot(1:min(Pat^2,disDim2),e(1:min(Pat^2,disDim2)));
    graph2d.constantline(0,'LineStyle',':','Color',[.7 .7 .7]);
    axis([1,min(Pat^2,disDim2),min(e),max(e)*0.9]);
    xlabel('Eigenvalue number');
    ylabel('Eigenvalue');


% Plot Correlation Coefficients
%     subplot(2,1,2);
%     hold on
%     corrs=[];
%     for i=1:5:min(size(Y1,2),Pat^2)
%         CORR = corrcoef(disMat,squareform(pdist(Y1(:,1:i))));
%         corrs=[corrs CORR(1,2)];
%     end
%     plot(1:5:min(size(Y1,2),Pat^2),corrs,'ko-','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerSize',8);
%     xlabel('Dimension');
%     ylabel('Correlation coefficient');


% Do MDS
    MDS = input('What is the MDS dimension ? ') ;
    maxerr = max(abs(pdist(X)-pdist(Y1(:,1:MDS)))) %#ok<NOPTS>
    Y = Y1(:,1:MDS);


% Plot final distances' correlation
%     MDStest = squareform(pdist(Y));
%     figure;
%     plot(disMat,MDStest,'b.');
%     CORR = corrcoef(disMat,MDStest);
%     xlabel('Dissimilarity distance');
%     ylabel(['Euclidean distance in ',num2str(MDS),' D']);
%     title(['Correlation Coefficient = ',num2str(CORR(1,2))]);

end
