function[avim] = averageImage(imar)
    nx = size(imar,1); ny = size(imar,2); nz = size(imar,3);
    sumim = zeros(nx,ny);
    for(i = 1:nz)
        sumim = sumim + imar(:,:,i);
    end
    avim = sumim ./ nz;
end