function[corr_coeff] = calculate_corr_vector(elems_to_compare,frstack)
    n = length(elems_to_compare);
    corr_coeff = zeros((n*(n-1))/2,2);
    %vectorize for speed
    ndx=1;
    for( i = 1:(n-2) )
        array1 = (frstack(:,:,i) - mean(mean(frstack(:,:,i))));
        [xx,yy] = arrayfun(@(ind) findmax(xcorr2_fft(sub_mean(frstack(:,:,i)),sub_mean(frstack(:,:,ind)) )),...
                                          (i+2):n);
      %  tic; [xx,yy] = arrayfun(@(ind) refine_autocorr_peak(xcorr2_fft(sub_mean(frstack(:,:,i)),sub_mean(frstack(:,:,ind)))),...
      %                             (i+1):n); toc;
        corr_coeff(ndx:(ndx+length(xx)-1),:) = [-(xx-size(frstack,1))' -(yy-size(frstack,2))'];
        ndx=ndx+length(xx);
    end
    corr_coeff = corr_coeff(1:(ndx-1),:);
end