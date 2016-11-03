function[rmsf] = atom_rmsf(pdb_ensemble)
    %assume pdbs are already aligned
    natom = length(pdb_ensemble{1}.X);
    rmsf = zeros(1,natom);
    npdb = length(pdb_ensemble);
    for(i = 1:natom)
        Xi = cellfun(@(thispdb)(thispdb.X(i)),pdb_ensemble);
        Yi = cellfun(@(thispdb)(thispdb.Y(i)),pdb_ensemble);
        Zi = cellfun(@(thispdb)(thispdb.Z(i)),pdb_ensemble);
        meanXi = mean(Xi);
        meanYi = mean(Yi);
        meanZi = mean(Zi);
        
        rmsf(i) = sqrt(sum([ (Xi-meanXi).^2 ,...
                             (Yi-meanYi).^2 ,...
                             (Zi-meanZi).^2 ])/npdb);
    end
end