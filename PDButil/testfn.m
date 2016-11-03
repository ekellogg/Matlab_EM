lst = textread('lst','%s','delimiter','\n');
pdbens = {};
for(i = 1:length(lst))
    pdbens{i} = readPDBnoH(lst{i});
    pdbens{i}.outfile = strcat(pdbens{i}.outfile,'flipd.pdb');
end

find_degenerate_aromaticpos(pdbens)
for(i = 1:length(pdbens))
    mat2pdb(pdbens{i});
end