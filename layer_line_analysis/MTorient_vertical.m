function[bestrot] = MTorient_vertical(im)
    %max = 0;
    rot = 0:1:179;
    proj_range = arrayfun(@(r_deg) (range(sum(imrotate(im,r_deg),1)))  ,rot);
    bestrot = rot( find(proj_range == max(proj_range)) );
   
    %vectorized code above simply does what is commented out below
    %for(i = 0:1:180)
    %   r = range(sum(imrotate(im,i),1));
    %   if( r  > max )
    %       max = r;
    %       rot = i;
    %   end
    %end
end