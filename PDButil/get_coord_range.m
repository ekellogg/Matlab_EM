function[range_r] = get_coord_range(inpdb);
    range_r = [min(inpdb.X) max(inpdb.X);
               min(inpdb.Y) max(inpdb.Y);
               min(inpdb.Z) max(inpdb.Z)];
end