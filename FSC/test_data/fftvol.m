function[F] = fftvol(vol)
    %F = fftshift( fftn( im ));
    F = fftshift( abs( fftn( im )));
    
end