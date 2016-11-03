function[res] = switch_atom_coordinates(res,at1,at2)
    [aa,ndx1] = find_atom_name(res,at1);
    [aa,ndx2] = find_atom_name(res,at2);
    tmp = [res.X(ndx1) res.Y(ndx1) res.Z(ndx1)];
    res.X(ndx1) = res.X(ndx2);
    res.Y(ndx1) = res.Y(ndx2);
    res.Z(ndx1) = res.Z(ndx2);
    res.X(ndx2) = tmp(1);
    res.Y(ndx2) = tmp(2);
    res.Z(ndx2) = tmp(3);
end