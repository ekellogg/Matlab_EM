function[binned_img] = bin_images(img,bin)
    binned_img = imresize(img,1/bin);
end