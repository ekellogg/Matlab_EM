function[summed_im] = apply_shifts_to_stack(shifts,frstack)
    summed_im = frstack(:,:,1);
%try resizing image to prevent interpolation artifacts
scale = 2;
      for(i = 1:size(shifts,1))
%         tic; summed_im = summed_im + gather(imresize( imtranslate( imresize( gpuArray(frstack(:,:,(i+1))), scale), shifts(i,:) * scale ), 1/scale)); toc

          tic; summed_im = summed_im + gather(imresize( imtranslate( imresize( gpuArray(frstack(:,:,1)), scale), shifts(i,:) * scale ), 1/scale)); toc
      end
end