function[atm,ndc] = find_atom_name(pdb,aname)
    ndc = find( cellfun(@(x)(strcmp(x,aname) == 1),pdb.atomName) );
    atm = get_ndx(pdb,ndc);
end