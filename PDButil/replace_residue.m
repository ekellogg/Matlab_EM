function[pdb] = replace_residue(pdb,res,pos)
    ndx = find(pdb.resNum == pos);
    pdb = replace_ndx(pdb,res,ndx);
end