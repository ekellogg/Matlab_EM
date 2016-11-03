addpath('~/projects/code/Matlab_EM/EMIODist/')
addpath('~/projects/code/Matlab_EM/MRCIO/')
addpath('~/projects/code/Matlab_EM/basic/')
addpath('code/');

fid = fopen('img_lst');
C= textscan(fid,'%s\n');
C=C{1};
apix=1.32;
outputdir = 'images/test_liz_simple_driftcorr_4';

for(i = 1:length(C))
    tic;
    mrc_i = squarify(readMRCfile(C{i}));
    fn = remove_path(C{i});
    
    dfcorr = simple_drift_correction(mrc_i);
    dfsum = cumsum(dfcorr,1);
    dfsum = gather(dfsum);
    dlmwrite(strcat(outputdir,'/',prefix(fn),'-drift.txt'),gather(dfcorr));
    summed_im = apply_shifts_to_stack_v2(dfsum,mrc_i,apix);
    
    WriteMRC(summed_im./size(mrc_i,3),apix,strcat(outputdir,'/',prefix(fn),'-simpledfcorr.mrc')); toc
end