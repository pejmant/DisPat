function boolOutput = existNonFrozenNodes(frozenRealiz, wx, wy, wz)

boolOutput = false;

frozenTemplate = frozenRealiz(wx, wy, wz);
if any(frozenTemplate(:) == 0)
    boolOutput = true;
end

end
