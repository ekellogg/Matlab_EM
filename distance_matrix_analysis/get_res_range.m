function[r] = get_res_range(index_start,index_end,pdb)
minn = min(pdb.resNum(index_start:index_end));
maxn = max(pdb.resNum(index_start:index_end));
r = [minn maxn];
end