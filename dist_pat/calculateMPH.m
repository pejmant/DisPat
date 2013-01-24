function hist1 = calculateMPH(realization)
% multiple point histogram of realization variable


    xmph = 4;
    zmph = 1;
    par.Dim = size(realization,1);
    par.Dimz = size(realization,3);
    disDim  = par.Dim-xmph +1;
    disDimz = par.Dimz-zmph +1;
    hist1 = zeros(1,1+bin2dec(num2str(ones(1,xmph^2*zmph))));


    for i=1:disDim
        for j=1:disDim
            for k=1:disDimz
                temp = 1+bin2dec(num2str(reshape(realization(i:i+xmph-1,j:j+xmph-1,k:k+zmph-1),1,[])));
                hist1(temp) = hist1(temp) + 1;
            end
        end
    end
%     figure1 = figure;
%     axes('Parent',figure1,'YScale','log','YMinorTick','on'); box('on'); hold('all');
%     semilogy(hist,'LineWidth',2,'Color',[0 0 0]);
%     xlabel('Multiple-point Configuration');
%     ylabel('Number of repetitions');
%     title('Proposed Method Realization');

end