function[alignedim] = alignImages(imset,params,mask)
    nx = size(imset,1); ny = size(imset,2); nz = size(imset,3);
    msize = maskSize(mask);
    mx = msize(1); my = msize(2);
    aimsz = [nx,ny];
    if( mx < nx ) 
       aimsz = [mx,my];  
    end
    alignedim = zeros(aimsz(1),aimsz(2));
    for(i = 1:nz)
       %frealign conventions:
       %particle origin in lower left corner
       % in-plane rotation PSI, positive = counter-clockwise rotation
       %to align with reference, a clockwise rotation (negative) needs to
       %be applied (same as imrotate, but imrotate rotates around image center)
       % X,Y shifts: positive values describe shifts right and up. negative
       % values need to be applied to align with reference
       % frealign applied shifts first and then rotation. no mirrors used
       % must apply upward shift because imgs are mirrored in matlab wrt
       % frealign
       
       %im = imrotate(imset(:,:,i),params(i,2));
       
       %OLD 
       im = imrotate(imshift(imset(:,:,i),-params(i,6),-params(i,5)),params(i,2));
       alignedim = alignedim + cropImage(maskImage(im,mask));
       if( mod(i,100) == 0 )
          display(sprintf('done with index %d',i)) 
       end
    end
    alignedim = alignedim ./ nz;
end