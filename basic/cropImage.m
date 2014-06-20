function[cim] = cropImage(im)
    [i,j] = find(im);
    cim = im(min(i):max(i),min(j):max(j));
end