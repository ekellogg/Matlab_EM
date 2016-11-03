function[summed_im] = apply_shifts_to_stack(shifts,frstack)
    summed_im = frstack(:,:,1);
%     for(i = 1:size(shifts,1))
%         tform = affine2d(make_translation_matrix(shifts(i,:)));
%         summed_im = summed_im + imwarp(frstack(:,:,1),tform,'linear');
%     end

%simplest way
%     for(i = 1:size(shifts,1))
%          summed_im = summed_im + imtranslate(frstack(:,:,(i+1)),shifts(i,:));
%     end

%try resizing image to prevent interpolation artifacts
scale = 4;
      for(i = 1:size(shifts,1))
         summed_im = summed_im + imresize( imtranslate( imresize( frstack(:,:,(i+1)), scale), shifts(i,:) * scale ), 1/scale);
      end
end