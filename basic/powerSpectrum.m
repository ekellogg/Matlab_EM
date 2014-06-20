function[ps] = powerSpectrum(stack)
    %expectation of abs-square fourier-transform
    nx = size(stack,1); ny = size(stack,2); nz = size(stack,3);
    avg = zeros(nx,ny);
    for(i = 1:nz)
        abs_fft = abs(fft(stack(:,:,i)));
        avg = avg + abs_fft*abs_fft;
    end
    ps = avg ./ nz;
end