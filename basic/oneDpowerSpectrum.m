function[freq,rotavg,r_ind] = oneDpowerSpectrum(im,apix)
    %rotationally averaged
    %adapted from appion/ace/getprofile_ang.m  
    imsize = size(im,1); %assumes square
    c = floor(size(im)/2)+1; %assumes center in middle of image
    absimfft = abs(fftim(im));
    %rotational averaging
    r_ind = [0:(round(imsize/2)-1)];
    theta_ind = 0:(2*pi/(6*c(1))):(2*pi);
    [r theta] = meshgrid(r_ind,theta_ind);
    %convert to cartesian coordinates
    i = c(1) + (r).*cos(theta);
    j = c(2) + (r).* sin(theta);
    i = i';
    j = j';
    val = interp2(absimfft,i,j,'bicubic');
    rotavg = mean(abs(val),2);
    freq = r_ind/(imsize*apix);
    r_ind = r_ind + c(1);
end