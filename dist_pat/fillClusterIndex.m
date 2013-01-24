function clusterIdx = fillClusterIndex(clusterIdx, idx)

for i = 1:max(idx)
    clusterIdx{i} = find(idx == i);
end
