function[nm] = remove_path(fn)
    alldirs = strfind(fn,'/');
    nm = fn((max(alldirs)+1):end);
end