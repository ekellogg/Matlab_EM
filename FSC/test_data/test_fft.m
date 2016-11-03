addpath('~/projects/code/Matlab_EM/basic/')
addpath('~/projects/code/Matlab_EM/MRCIO/')
addpath('~/projects/code/Matlab_EM/EMIODist/')

map1 = ReadMRC('14pf_r0_map1.mrc');
map2 = ReadMRC('14pf_r0_map2.mrc');

F = fsc(map1,map2);