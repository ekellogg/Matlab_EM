function[pdb_ens,rmsf,ndx] = replace_bfact(pdblst,nstd_dev)
    lst = textread(pdblst,'%s','delimiter','\n');
    pdb_ens = {};
    for(i = 1:length(lst))
        pdb_ens{i} = readPDBnoH(lst{i});
    end
    %get rid of aromatic flips.. they are equivalent!
    pdb_ens = find_degenerate_aromaticpos(pdb_ens);
    rmsf = atom_rmsf(pdb_ens);
    ndx = find( rmsf < (mean(rmsf)+nstd_dev*std(rmsf)) );
    %cap everything above threshold at max b-factor,
    %scale everything else linearly from the lowest value
end