function[uniqchains] = all_chains(pdb)
    uniqchains = unique(pdb.chainID);
end