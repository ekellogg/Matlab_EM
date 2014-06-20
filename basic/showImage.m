function[handl] = showImage(im)
    minVal = min(min(im));
    im = im - minVal;
    maxVal = max(max(im));
    image( im / maxVal * 100 );
    
end