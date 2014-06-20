calim = zeros(7676,7420,3);
calim(:,:,1) = ReadMRC('images/calibration_images_061314/SR_focus.mrc');
calim(:,:,2) = ReadMRC('images/calibration_images_061314/SR_1p75.mrc');
calim(:,:,3) = ReadMRC('images/calibration_images_061314/SR_3.mrc');

[fr1,pw1d1,r_ind1] = oneDpowerSpectrum(calim(1:7420,:,1),1);
smo1 = smoothn(pw1d1,31);
[fr2,pw1d2,r_ind2] = oneDpowerSpectrum(calim(1:7420,:,2),1);
smo2 = smoothn(pw1d2,31);
[fr3,pw1d3,r_ind3] = oneDpowerSpectrum(calim(1:7420,:,3),1);
smo3 = smoothn(pw1d3,31);

gold_band_ind_begin = 5650;
gold_band_ind_end = 5950;

search_band_ind_begin = find(r_ind1 == 5650);
search_band_ind_end = find(r_ind1 == 5950);

search_reg = smo1(search_band_ind_begin:search_band_ind_end);

max_ind1 = r_ind1( search_band_ind_begin + find(search_reg == max(search_reg)) - 1 );

search_band_ind_begin = find(r_ind2 == 5650);
search_band_ind_end = find(r_ind2 == 5950);

search_reg = smo2(search_band_ind_begin:search_band_ind_end);

max_ind2 = r_ind2( search_band_ind_begin + find(search_reg == max(search_reg)) - 1 );

search_band_ind_begin = find(r_ind3 == 5650);
search_band_ind_end = find(r_ind3 == 5950);

search_reg = smo3(search_band_ind_begin:search_band_ind_end);

max_ind3 = r_ind3( search_band_ind_begin + find(search_reg == max(search_reg)) - 1 );

gold_line = 2.3469 % gold lines show up at 2.3 Angstrom
im_size = 7420;
im_center = 3710; %7420/2

apix1 = 1/((1/gold_line)*(im_size)/(max_ind1-im_center))*2; %super-resolution mode is 2x resolution of counting mode
apix2 = 1/((1/gold_line)*(im_size)/(max_ind2-im_center))*2;
apix3 = 1/((1/gold_line)*(im_size)/(max_ind3-im_center))*2;

display(sprintf('ang/pix for image 1: %d, image 2: %d, image 3: %d, mean: %d',apix1,apix2,apix3,mean([apix1,apix2,apix3])))