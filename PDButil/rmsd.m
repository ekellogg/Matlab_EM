function[r] = rmsd(pdb1,pdb2)
    r = sqrt(sum((pdb1.X-pdb2.X).^2 + (pdb1.Y-pdb2.Y).^2 + (pdb1.Z-pdb2.Z).^2)/length(pdb1.X));
end