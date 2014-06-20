function[sz] = maskSize(m)
    [i,j] = find(m);
    sz_x = max(i) - min(i) + 1;
    sz_y = max(j) - min(j) + 1;
    sz = [sz_x sz_y];
end