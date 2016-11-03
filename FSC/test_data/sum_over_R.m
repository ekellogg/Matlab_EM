function[s] = sum_over_R(f)
    imsize = length(f);
    c = floor(length(imsize/2))+1;
    r_ind = [0:(round(imsize/2)-1)];
    theta_ind = 0:(2*pi/100):(2*pi);
    [r theta] = meshgrid(r_ind,theta_ind);
    i = c + (r).*cos(theta);
    j = c + (r).*sin(theta);
    v = interp2(f,i,j,'bicubic');
    s = sum(v,1);
%    freq = r_ind/(imsize*apix);
end