addpath('~/projects/code/Matlab_EM/EMIODist/')
addpath('~/projects/code/Matlab_EM/MRCIO/')
addpath('~/projects/code/Matlab_EM/basic/')
addpath('drift_correction_code/');

%load example image
p = imread('pout.tif');
sp = squarify(p);

%generate fake shifts
%fake_shifts = [0:1:10; (0:1:10).^2]';
fake_shifts = fake_shifts(1:4,:)

%apply fake shifts to image to create fake trajectory
frstack = zeros([size(sp) size(fake_shifts,1)+1]);
frstack(:,:,1) = sp;
for(i = 1:size(fake_shifts,1))
    frstack(:,:,i+1) = imshift_fft((fake_shifts(i,:)),sp);
end

%apply drift correction to get back the shifts
dfcorr = simple_drift_correction(frstack);
dfsum = gather(cumsum(dfcorr,1));

clf
%invert the shifts to get back a full image
subplot(2,2,1)
plot(fake_shifts(:,1),fake_shifts(:,2));
hold on
plot([0 dfsum(:,1)'],[ 0 dfsum(:,2)'],'-red');
legend('input','output')

%this function will invert shifts so that you don't have to.
%dfsum(:,2) = -1*dfsum(:,2);
dfsumflipd = [dfsum(:,2) dfsum(:,1)];
summed_im = apply_shifts_to_stack_v2(dfsumflipd,frstack);
subplot(2,2,2)
showImage(summed_im);
subplot(2,2,3);
showImage(sp)
subplot(2,2,4);
showImage(apply_shifts_to_stack_v2(dfsumflipd,frstack))
% showImage(abs(summed_im-sp))

