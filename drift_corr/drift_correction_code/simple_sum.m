function[avgim] = simple_sum(frames)
avgim = frames(:,:,1);
if( size( frames,3 ) > 1 )
    for(i = 2:size(frames,3))
        avgim = avgim + frames(:,:,i);
    end
    avgim = avgim ./ size(frames,3);
end
end