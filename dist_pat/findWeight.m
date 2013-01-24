function wieghtEvent = findWeight(hardDataMoved, frozenRealiz, wx, wy, wz)
% gives weight to the hard data (0.5), inner patch (0.3) and 0.2 to the rest

hardDataMoved = hardDataMoved(wx, wy, wz);
frozenRealiz  = frozenRealiz (wx, wy, wz);


wieghtEvent = 0.2*ones(size(frozenRealiz));
wieghtEvent(frozenRealiz ==1     ) = 0.3;
wieghtEvent(~isnan(hardDataMoved)) = 0.5;

wieghtEvent = reshape(wieghtEvent,1,[]);

end