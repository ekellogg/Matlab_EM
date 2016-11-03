function [chainout ] = get_chain(pdbin,chain)
    chainout = get_ndx(pdbin,find_chain(pdbin,chain));
end