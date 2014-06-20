function mc=Crop(m,n,isstack,fillval)
% function mc=Crop(m,n,isstack,fillvalue)
% Reduce the size of the 1d array, square or cube m by cropping
% (or increase the size by padding with fillval, by default zero)
% to a final size of n x n or n x n x n.  This is the analogue of
% Downsample, but doesn't change magnification.
% This handles odd and even-sized arrays correctly  The center of an
% odd array is taken to be at (n+1)/2, and an even array is n/2+1.
% If the flag isstack = 1 then a 3D array m is treated as a stack of 2D
% images, and each image is cropped to n x n.
% For 2D images, the input image doesn't have to be square.

if nargin<3
    isstack=0;
end;
if nargin<4
    fillval=0;
end;

sz=size(m);
ndi=ndims(m);
if ndi==2 && any(sz==1)
    ndi=1;
end;

switch ndi
    case 1
        n1=numel(m);
        m=reshape(m,n1,1); % force a column vector
        ns=floor(n1/2)-floor(n/2);  % Shift term for scaling down.
        if ns>=0 % cropping down
            mc=m(ns+1:ns+n);
        else
            mc=fillval*ones(n,1);
            mc(1-ns:n1-ns)=m;
        end;
              
    case 2
        nx=size(m,1);
        ny=size(m,2);
        nsx=floor(nx/2)-floor(n/2);  % Shift term for scaling down.
        nsy=floor(ny/2)-floor(n/2);
        if nsx>=0 % cropping down
            mc=m(nsx+1:nsx+n,nsy+1:nsy+n);
        else
            mc=fillval*ones(n,n);
            mc(1-nsx:nx-nsx,1-nsy:ny-nsy)=m;
        end;
        
    case 3 % m is 3D
        [n1 n2 ni]=size(m);
        ns=floor(n1/2)-floor(n/2);  % Shift term for scaling down.
        
        if isstack % a stack of 2D images
            if ns>=0 % cropping down
                mc=m(ns+1:ns+n,ns+1:ns+n,:);
            else
                mc=fillval*ones(n,n,ni);
                mc(1-ns:n1-ns,1-ns:n1-ns,:)=m;
            end;
        else
            if ns>=0 % cropping down
                mc=m(ns+1:ns+n,ns+1:ns+n,ns+1:ns+n);
            else
                mc=fillval*ones(n,n,n);
                mc(1-ns:n1-ns,1-ns:n1-ns,1-ns:n1-ns)=m;
            end;
        end;
    otherwise
        error(['Downsample: dimension too large: ' num2str(ndi)]);
end;
