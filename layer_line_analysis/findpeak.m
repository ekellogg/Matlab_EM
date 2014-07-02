function[xval] = findpeak(x,y,min_x,max_x)
    xend = min(find(x < min_x));
    xstart = max(find( x > max_x ));
    y_seg = y(xstart:xend);
    x_seg = x(xstart:xend);
    max_ndx = find(y_seg == max(y_seg));
    xval = mean(x_seg(max_ndx));
end