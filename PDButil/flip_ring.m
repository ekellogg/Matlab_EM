function[flipdres] = flip_ring(res)
    flipdres = res;
    flipdres = switch_atom_coordinates(flipdres,'CE1','CE2');
    flipdres = switch_atom_coordinates(flipdres,'CD1','CD2');
end