function [dataEvent, status] = getDataEvent(realization, wx, wy, wz)

% gets the data event located at location wx and wy
% also checks if the dataevent is empty (all NaNs) or it has
% some data conditioning
Pat  = length(wx);
Patz = length(wz);

dataEvent = reshape(realization(wx,wy,wz), 1, []);

% DT=bwdist(realization(wx,wy));    % not for 3D case
% if DT==Inf
%     DT = ones(1,Pat^2);
% end
% dataEvent=reshape(1-DT./max(max(DT)),1,Pat^2);



% check HOW many grid nodes are unvisited
diff = sum(abs(dataEvent - 0.5));

switch diff
    case 0
        status = 'empty';
    case Pat^2*Patz
        status = 'full';
    otherwise
        status = 'some';
end


end