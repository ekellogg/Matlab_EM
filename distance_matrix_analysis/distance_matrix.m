function[dmat] = distance_matrix(pdb)
    nat = length(pdb.X);
    dmat = zeros(nat,nat);
    for(i = 1:nat)
        for(j = 1:(i-1))
            dd = sqrt((pdb.X(i) - pdb.X(j)).^2 + (pdb.Y(i)-pdb.Y(j)).^2 + (pdb.Z(i)-pdb.Z(j)).^2) ;
            dmat(i,j) = dd;
            dmat(j,i) = dd;
        end
    end
end