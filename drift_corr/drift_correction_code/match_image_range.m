function[newim] = match_image_range(imtgt,imref)
    [rmin_ref,rmax_ref] = image_range(imref);
    [rmin_tgt,rmax_tgt] = image_range(imtgt);
    
    norm_imtgt = normalize_image(imtgt);
    norm_imref = normalize_image(imref);
    
    newim = (norm_imtgt-rmin_tgt) ./ (rmax_tgt-rmin_tgt) * (rmax_ref - rmin_ref) + rmin_ref;
    %newim = (imtgt - rmin_tgt) ./ (rmax_tgt - rmin_tgt) * rmax_ref;
end