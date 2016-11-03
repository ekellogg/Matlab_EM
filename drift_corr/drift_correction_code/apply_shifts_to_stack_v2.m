function[summed_im] = apply_shifts_to_stack(shifts,frstack)

    %do this using ffts instead
    summed_im = gpuArray(fftshift(fft2(frstack(:,:,1))));
    sz = size(frstack);
    %define frequency across fourier domain
    xs = -sz(1)/2; xe = (sz(1)/2)-1;
    ys = -sz(2)/2; ye = (sz(2)/2)-1;
    [xf,yf] = meshgrid( (xs:1:xe)./sz(1) , (ys:1:ye)./sz(2) ); %mesh of fourier frequencies
    xf = gpuArray(xf);
    yf = gpuArray(yf);
    
    for(i = 1:size(shifts,1))
       summed_im = summed_im + fftshift( fft2( gpuArray(frstack(:,:,(i))) )) ...
                               .*exp(-1i*2*pi.*(xf*-shifts(i,1) + yf*-shifts(i,2)) ); 
    end
    summed_im = gather(real(ifft2(ifftshift( summed_im  ))));
end