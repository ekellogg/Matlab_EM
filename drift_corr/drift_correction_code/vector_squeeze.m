function[v] = vector_squeeze(im)
    v = reshape(im,1,size(im,1)*size(im,2));
end