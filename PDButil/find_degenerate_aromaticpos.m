function[] = find_degenerate_aromaticpos(pdbens)
 %find all tyrosines and phenylalanines
 %minimize the rmsd to the first pdb
    resname = pdbens{1}.resName;
    alltyr = find(cellfun('isempty',cellfun(@(x)(strfind(x,'TYR')),resname,'UniformOutput',0)) == 0);
    resnum = unique(pdbens{1}.resNum(alltyr));
    for(i = 1:length(resnum))
        ref_pdb_i = pdbens{1};
        for(j = 2:length(pdbens))
            
        end
    end
    
    allphe = find(cellfun('isempty',cellfun(@(x)(strfind(x,'PHE')),resname,'UniformOutput',0)) == 0);
    
end