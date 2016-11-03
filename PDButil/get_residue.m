function[res] = get_residue(pdb,pos)
    %resturns a pdb like data structure but containing only a single
    %residue
    ndc = find(pdb.resNum == pos);
    res = get_ndx(pdb,ndc);
end