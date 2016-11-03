function[ndx] = find_chain(pdb,ch)
 ndx = find(cellfun(@(x) strcmp(x,ch) == 1,pdb.chainID,'UniformOutput',1 ));
end