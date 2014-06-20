function[mim] = maskImage(im,mask)
    imsize = size(im,1);
    masksize = size(mask,1);
    if( imsize > masksize )
        m = zeros(imsize,imsize);
        ori = ceil((imsize - masksize)/2);
        if( ori == 0 )
            pause(1)
        end
        m(ori:(ori+masksize-1),ori:(ori+masksize-1)) = mask;
        mask = m;
    end
        mim = im .* mask;
end