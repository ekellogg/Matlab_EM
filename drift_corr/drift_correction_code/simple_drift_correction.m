function[img_shift] = simple_drift_correction(fr)
    %input frames are likely not gain corrected.
    BFACT = 150; 
    for(i = 1:size(fr,3))
       fr(:,:,i) = apply_bfact(fr(:,:,i),BFACT); 
    end
    nfr = size(fr,3);
    %create matrix A
    A = create_coeff_matrix(1:nfr);
    b = calculate_corr_vector_v2(1:nfr,fr);
    img_shift = inv(transpose(A)*A)*transpose(A)*b;
    
    %r = inv(transpose(A)*A)*transpose(A)*b;
    %residual_error = sum(((A*r) - b).^2,2);
    %coeff_to_keep = find(residual_error < 1);
    %A = A(coeff_to_keep,:); 
    %img_shift = inv(transpose(A)*A)*transpose(A)*b(coeff_to_keep,:); %maybe this is incorrect?
end