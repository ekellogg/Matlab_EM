function[] = show3dimage(mrc)
subplot(1,3,1);
showImage(squeeze(sum(mrc,1)));
subplot(1,3,2);
showImage(squeeze(sum(mrc,2)));
subplot(1,3,3);
showImage(squeeze(sum(mrc,3)));
colormap(gray(100))
end
