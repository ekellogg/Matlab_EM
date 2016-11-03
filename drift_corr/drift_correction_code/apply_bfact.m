function[bimg] = apply_bfact(img,bfactor)
%    bimg = ifft2(fft2(img).*exp(-0.5*bfactor));
    %add padding in FFT to minimize interpolation error
    c = ceil(size(img))./2;
    x = 1:1:size(img,1);
    y = 1:1:size(img,2);
    [X,Y] = meshgrid(x,y);
%    f = sqrt((X-c(1)).^2 + (Y-c(2)).^2)./(size(img,1)*apix); %convert to frequency
    %what the fuck do you use here for spatial frequency?? 1/d is not
    %right because it ends up amplifying high resolution features, not low
    %resolution features. 
    %f = 1./(sqrt((X-c(1)).^2 + (Y-c(2)).^2));
    %above is not right
    f = (sqrt((X-c(1)).^2 + (Y-c(2)).^2))./(size(img,1));
    z = exp(-0.5*bfactor*(f.^2));
    bimg = real(ifft2(ifftshift(  fftshift(fft2(img)).*z  )));
end