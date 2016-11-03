function[A] = create_coeff_matrix(i)
    A = gpuArray(zeros(length(i),length(i)-1));
    ii = i(1);
    ndx = 1;
    while(ii <= i(length(i)))
        jj = ii+1;
        while(jj < i(length(i))  )
            A(ndx,ii:jj) = 1;
            %do not align frames that are adjacent
            jj = jj + 1;
            ndx = ndx + 1;
        end
        ii = ii + 1;
    end
end