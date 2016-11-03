function[corr_coeff] = calculate_corr_vector(elems_to_compare,frstack)
    n = length(elems_to_compare);
    corr_coeff = zeros((n*(n-1))/2,2);
    ndx = 1;
    for( i = 1:n )
        for( j = (i+1):(n-1) )
            %convert to gpu arrays
            array1 = (frstack(:,:,i) - mean(mean(frstack(:,:,i))));
            array2 = (frstack(:,:,j) - mean(mean(frstack(:,:,j))));
            tic; autocorr_map = (xcorr2_fft(array1,array2));
            sprintf('finished computing auto-corr between elems %d and %d',i,j); toc
            %best x and y
            %max index
%            autocorr_map(size(frstack,1),size(frstack,2)) = 0;
            [val,index] = max(autocorr_map(:));
            [x,y] = ind2sub(size(autocorr_map),index);
            x = x - size(frstack,1);
            y = y - size(frstack,2);
            corr_coeff(ndx,:) = [x,y];
            ndx = ndx + 1;
        end
    end
end