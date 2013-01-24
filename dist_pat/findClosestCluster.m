function idxNumber = findClosestCluster(dataEvent, prototype, wieghtEvent)

% This function finds the closest cluster to the dataevent

index = find(dataEvent ~= 0.5);

% only search through Non-NaN s
reducedPrototype = prototype(:, index);
reducedDataEvent = dataEvent(1, index);
wieghtEvent      = wieghtEvent(1,index);

difference = sum(bsxfun(@times,abs(bsxfun(@minus, reducedPrototype, reducedDataEvent)),wieghtEvent), 2);
[dummy, idxNumber] = min(difference);

end
