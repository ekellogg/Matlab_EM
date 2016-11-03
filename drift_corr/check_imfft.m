% s = imread('pout.tif');
% sp = squarify(s);
% clf
% subplot(1,2,1)
% tt = imshift_fft([6 6],sp,1);
% tt = tt - min(tt(:)); %normalize image range
% image( tt  ./ max(tt(:)) * 100 ); 
% subplot(1,2,2)
% tt = imshift_fft([32 32],sp,1);
% tt = tt - min(tt(:));
% image( tt ./ max(tt(:)) * 100);
clf
fake_shifts = [(0:1:50); (0:1:50)]';

cols = jet(length(fake_shifts));
for(i = 1:size(fake_shifts,1))
subplot(1,2,1)
showImage(imshift_fft(fake_shifts(i,:),sp,1))
title(sprintf('%d pi',fake_shifts(i,:)))
subplot(1,2,2)
hold on
plot(sort(reshape(imshift_fft(fake_shifts(i,:),sp,1),[ 240*240 1])),'Color',cols(i,:))
pause(0.5)
end