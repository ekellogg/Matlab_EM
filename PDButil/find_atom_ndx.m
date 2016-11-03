function[ndx] = find_atom_ndx(pdb,chainid,resnum,atomname)
    chain_ii = find(cellfun(@(x) strcmp(x,chainid) == 1,pdb.chainID,'UniformOutput',1 ));
    resnum_ii = find( pdb.resNum == resnum);
    atomname_ii = find(cellfun(@(x) strcmp(x,atomname) == 1,pdb.atomName,'UniformOutput',1));
    ndx = intersect(intersect(chain_ii,resnum_ii),atomname_ii);
end