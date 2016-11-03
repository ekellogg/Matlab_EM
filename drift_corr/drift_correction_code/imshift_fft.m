function[shiftim] = imshift_fft(shift,fr)
sz = size(fr);
%mesh of fourier frequencies
[xf,yf] = meshgrid( ((-sz(1)/2):1:(sz(1)/2)-1)./(sz(1)) , ((-sz(2)/2):1:(sz(2)/2)-1)./(sz(2)) ); 
F = fftshift( fft2( fr ));

%define frequency throughout fourier space
%Fshift=  F.*exp(-1i*2*pi.*(xf*shift(1) + yf*shift(2)) / (sz(1)) );

%changed according to page 322 in frank's EM book.
Fshift=  F.*exp(-1i*(2*pi.*(xf*shift(1) + yf*shift(2))) );
shiftim = real(ifft2(ifftshift( Fshift  )));

%shiftim = match_image_range(shiftim,fr);
end