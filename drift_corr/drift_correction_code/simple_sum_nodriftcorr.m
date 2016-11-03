addpath('~/projects/code/Matlab_EM/EMIODist/')
addpath('~/projects/code/Matlab_EM/MRCIO/')
addpath('~/projects/code/Matlab_EM/basic/')
addpath('code/');

fid = fopen('img_lst');
C= textscan(fid,'%s\n');
C=C{1};
apix=1.32;
outputdir = 'images/control_nodrift_correction';

for(i = 1:length(C))
    tic;
    mrc_i = squarify(readMRCfile(C{i}));
    fn = remove_path(C{i});
    avg_i = sum(mrc_i,3)./size(mrc_i,3);
    WriteMRC(avg_i,apix,strcat(outputdir,'/',prefix(fn),'-avg.mrc')); toc
end
