function[newim] = fft_crop(im,scale)
    newsize = size(im)/scale;
    c = size(im)./2;
    s = c - (newsize(1)/2) + 1;
    e = c + (newsize(2)/2);
    fftim = fftshift(fft2(im));
    newim = real(ifft2(ifftshift( fftim( s(1):e(1), s(2):e(2) ))));
end