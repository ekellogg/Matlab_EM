function[x,ll] = get_layer_lines(im,apix)
    fft_im = fftim(im);
    ll = sum(fft_im,1);
    nx = size(im,1);
    ori = (nx/2)+1; %indexed from 1
    x = nx*apix./((ori:length(ll)) - ori);
    ll = ll(ori:end);
end