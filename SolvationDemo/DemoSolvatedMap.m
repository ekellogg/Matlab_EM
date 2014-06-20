pdbPath='/home/ekellogg/projects/taxol_MT/simulated_maps/';
[coords types]=ReadPDBAtoms([pdbPath 'extended_K.pdb']);
pixA=1.32;  % angstroms per pixel

[TotalDens volm]=SolventAndProteinDensity(coords, types);
n1=size(TotalDens,1);
n0=2*ceil(n1/(2*pixA));  % size of final map, even to avoid possibe bugs

%% Resample simulations to the data pixA value
TotalDens=TotalDens-TotalDens(1);  % remove DC
sim=DownsampleGeneral(TotalDens,n0,1/pixA)*pixA; % Convert from 1A/pixel
ShowSections(sim);
