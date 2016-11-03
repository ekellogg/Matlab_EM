function[res] = get_ndx(pdb,ndx)
    res = new_res();
    res.outfile = pdb.outfile;
    res.recordName = {pdb.recordName{ndx}};
    res.atomNum = pdb.atomNum(ndx);
    res.atomName = {pdb.atomName{ndx}};
    res.altLoc = {pdb.altLoc{ndx}};
    res.resName = {pdb.resName{ndx}};
    res.chainID = {pdb.chainID{ndx}};
    res.resNum = pdb.resNum(ndx);
    res.X = pdb.X(ndx);
    res.Y = pdb.Y(ndx);
    res.Z = pdb.Z(ndx);
    res.occupancy = pdb.occupancy(ndx);
    res.betaFactor = pdb.betaFactor(ndx);
    res.element = {pdb.element{ndx}};
    res.charge = {pdb.charge{ndx}};
end