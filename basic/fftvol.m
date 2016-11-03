function[F] = fftvol(vol)
    F = fftshift( abs( fftn( vol )));
    %F = fftshift(fft2(im)); %??
end