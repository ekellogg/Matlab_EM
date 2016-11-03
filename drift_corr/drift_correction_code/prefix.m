function[pr] = prefix(fn)
    pos = strfind(fn,'.');
    pr = fn(1:(pos-1));
end