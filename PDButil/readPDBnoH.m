function[pdb] = readPDBnoH(fn)
    pdb = pdb2mat(fn);
    %filter out hydrogens
    aname = pdb.atomName;
    Hs = cellfun(@(x)strfind(x,'H'),aname,'UniformOutput',0);
    heavyelems=find(cellfun('isempty',Hs));
    pdb.recordName = {pdb.recordName{heavyelems}};
    pdb.atomNum = [pdb.atomNum(heavyelems)];
    pdb.atomName = {pdb.atomName{heavyelems}};
    pdb.altLoc = {pdb.altLoc{heavyelems}};
    pdb.resName = {pdb.resName{heavyelems}};
    pdb.chainID = {pdb.chainID{heavyelems}};
    pdb.resNum = [pdb.resNum(heavyelems)];
    pdb.X = [pdb.X(heavyelems)];
    pdb.Y = [pdb.Y(heavyelems)];
    pdb.Z = [pdb.Z(heavyelems)];
    pdb.occupancy = [pdb.occupancy(heavyelems)];
    pdb.betaFactor = [pdb.betaFactor(heavyelems)];
    pdb.element = {pdb.element{heavyelems}};
    pdb.charge = {pdb.charge{heavyelems}};
end