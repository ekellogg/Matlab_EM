function[xx,yy] = refine_autocorr_peak(autocorr_map)
    [max_val,max_ind] = max(autocorr_map(:));
    [x,y]=ind2sub(size(autocorr_map),max_ind);
    ssize = 16;
    smap = autocorr_map( (x-(ssize/2)+1):(x+(ssize/2)), ...
                         (y-(ssize/2)+1):(y+(ssize/2)) );
    scaled_map = pad_fft(smap,512/16);
    scaledsize = size(scaled_map);
    origsize = size(autocorr_map);
    
    [smax,sind] = max(scaled_map(:));
    [sx,sy] = ind2sub(size(scaled_map),sind);
    xx = (sx - (scaledsize(1)/2))*(ssize/scaledsize(1)) + x;
    yy = (sy - (scaledsize(2)/2))*(ssize/scaledsize(2)) + y;
end