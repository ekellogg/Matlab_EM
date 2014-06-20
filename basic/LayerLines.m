function[x,ll] = LayerLines(im,apix)
    %assumes oriented vertically:
    %   |       | 
    %   |       |
    %   |       |
    ll=sum(fftim(im),2);
    origin = ceil(size(im,1)/2);
    ll=ll(origin:-1:1);
    x=(size(im,1)*apix)./(1:origin);
end