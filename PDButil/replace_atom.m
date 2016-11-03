function[newres] = replace_atom(res,aname,newatm)
    [at,nn] = find_atom_name(res,aname);
    newres = replace_ndx(res,newatm,nn);
end