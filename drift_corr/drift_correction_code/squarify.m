function [ sqfr ] = squarify( fr )
mindim = min(size(fr,1),size(fr,2));
sqfr = zeros(mindim,mindim,size(fr,3));
for(i = 1:size(fr,3))
    sqfr(:,:,i) = fr(1:mindim,1:mindim,i);
end

end

