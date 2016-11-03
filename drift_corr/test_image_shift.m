addpath('~/projects/code/Matlab_EM/EMIODist/')
addpath('~/projects/code/Matlab_EM/MRCIO/')
addpath('~/projects/code/Matlab_EM/basic/')
addpath('drift_correction_code/');
%load example image
p = imread('pout.tif');
sp = squarify(p);
fake_shifts = fake_shifts(1:4,:);