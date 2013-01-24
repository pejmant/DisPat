function hardDataMoved = moveHardData(hardData, par)
% move the hard data to nearest location

hardDataMoved = NaN(size(hardData));

wx = 1:par.m1:size(hardData,1);
wy = 1:par.m1:size(hardData,2);
wz = 1:par.m1:size(hardData,3);

[h1,h2,h3] = find(~isnan(hardData));

for i = 1:length(h1)
    [minh1, indminh1] = min(abs(bsxfun(@minus,wx,h1(i))));
    [minh2, indminh2] = min(abs(bsxfun(@minus,wy,h2(i))));
    [minh3, indminh3] = min(abs(bsxfun(@minus,wz,h3(i))));
    
    hardDataMoved(wx(indminh1(1)),wy(indminh2(1)),wz(indminh3(1))) = hardData(h1(i),h2(i),h3(i));
end


end