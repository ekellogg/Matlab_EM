function[res] = get_res(pdb,resnum)
    res.outfile = pdb.outfile;
    res_ndx = find(pdb.resNum == resnum);
    res.recordName = pdb.recordName{res_ndx};
    res.atomNum = pdb.atomNum{res_ndx};
    res.altLoc = pdb.altLoc{res_ndx};
    res.resName = pdb.resName{res_ndx};
    
end