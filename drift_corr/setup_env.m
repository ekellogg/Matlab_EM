addpath('~/projects/code/Matlab_EM/EMIODist/')
addpath('~/projects/code/Matlab_EM/MRCIO/')
addpath('~/projects/code/Matlab_EM/basic/')
addpath('drift_correction_code/');

 fr = readMRCfile('test_image/14nov18a_2_00041gr_00037sq_v01_00003hl_00002en.frames.mrc');
 fr = squarify(fr);
 lowpassfr = zeros(size(fr));
 for(i = 1:size(fr,3))
     lowpassfr(:,:,i) = medfilt2(fr(:,:,i),[5 5]);
 end

cropped_fr = zeros(ceil(size(fr)./[2 2 1]));
center_fr = ceil([size(fr,1)/2 size(fr,2)/2]);
for(i = 1:size(fr,3))
   cropped_fr(:,:,i) = fr( (center_fr(1)+1):end,...
                           (center_fr(2)+1):end,...
                           i); 
end

croppedfilt = zeros(size(cropped_fr));
for(i = 1:size(cropped_fr,3))
   croppedfilt(:,:,i) = medfilt2(cropped_fr(:,:,i),[5 5]); 
end

testimage = zeros(100,100);
for( i = 25:75 )
   testimage(i,50) = 1;
   testimage(50,i) = 1;
end