function[outx,outy] = findmax(tmp)
    [mm,ii] = max(tmp(:));
    [outx,outy] = ind2sub(size(tmp),ii);
end