function[pkloc_80A,peak_loc,layer_lines] = measure_80A_layer_lines(imgs,apix,padfactor)
    nz = size(imgs,3);
    rot = zeros(1,nz);
    pkloc_80A = zeros(1,nz);
    for(i = 1:nz)
        %create an image 2x the size to create finer sampling in fourier
        %space
        rot(i) = MTorient_vertical(imgs(:,:,i));
    end
    sprintf('done computing optimal rotations')
    nx  = size(imgs,1);
    ny  = size(imgs,2);
    peak_search_start = 65;
    peak_search_end = 90;
    peak_loc = {};
    layer_lines = {};
    for(i = 1:nz)
        if( mod(i,100) == 0 )
            sprintf('done with particle %d',i)
        end
        impx = zeros(nx*padfactor,nx*padfactor);
        lb = (nx*padfactor/2)-(nx/2);
        ub = (nx*padfactor/2)+(nx/2)-1;
        impx(lb:ub,lb:ub) = imgs(:,:,i);
        [x,ll] = get_layer_lines(imrotate(impx,rot(i)+90,'bilinear','crop'),apix);
        startndx = min(find( x <= peak_search_start ));
        endndx = max(find( x >= peak_search_end));
        xseg = x(endndx:startndx);
        llseg = ll(endndx:startndx);
        
        peak_loc{i} = x;
        layer_lines{i} = ll;
        pkloc_80A(i) = xseg( find(llseg == max(llseg)) );
    end
end