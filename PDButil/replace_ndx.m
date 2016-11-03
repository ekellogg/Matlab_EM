function[pdb] = replace_ndx(pdb,res,ndx)
for(nn = 1:length(ndx))
    n = ndx(nn);
    pdb.recordName{n} = res.recordName{nn};
    pdb.atomNum(n) = res.atomNum(nn);
    pdb.atomName{n} = res.atomName{nn};
    pdb.altLoc{n} = res.altLoc{nn};
    pdb.resName{n} = res.resName{nn};
    pdb.chainID{n} = res.chainID{nn};
    pdb.resNum(n) = res.resNum(nn);
    pdb.X(n) = res.X(nn);
    pdb.Y(n) = res.Y(nn);
    pdb.Z(n) = res.Z(nn);
    pdb.occupancy(n) = res.occupancy(nn);
    pdb.betaFactor(n) = res.betaFactor(nn);
    pdb.element{n} = res.element{nn};
    pdb.charge{n} = res.charge{nn};
end
end