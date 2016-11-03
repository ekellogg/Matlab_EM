function[F] = fsc(map1,map2)
    fftmap1 = fftn(map1);
    fftmap2 = fftn(map2);
    %need to sum over all r, currently these are 2D nxn matrices
    %a = sum_over_R(squeeze(dot(fftmap1,fftmap2)));
    %b = sum_over_R(squeeze(dot(fftmap1,fftmap1)));
    %c = sum_over_R(squeeze(dot(fftmap2,fftmap2))); %after summing over R, should be 1D:
    % 1 x n arrays.
    a = squeeze(dot(fftmap1,fftmap2));
    b = squeeze(dot(fftmap1,fftmap1));
    c = squeeze(dot(fftmap2,fftmap2));
    
    F = a ./ sqrt(b.*c);
end