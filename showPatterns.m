% By clicking on points in graph, you can later on
% input the numbers to see the coresponding pattern

if MDS==2 || PC==2
    yesno = input('Do you want to see some patterns related to the plot ? (Yes=1 / No=0) ');
    if yesno==1
        hold on
        if MDS==2
            plot(Y(:,1),Y(:,2),'.w')
        else
            plot(scores(:,1),scores(:,2),'.w')
        end
        gname;
        hold off
        no = input('input the corresponding case number ? (-1 to quit) ');
        while (no ~= -1)
            figure;
            imagesc(reshape(X(no,:),Pat,Pat))
            no = input('input the corresponding case number ? (-1 to quit) ');
        end
    end
end
