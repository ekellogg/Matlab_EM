function[im] = readDMMRC(filename)
    [fid,msg]=fopen(filename,'r');
    status=fseek(fid,0,-1);
    a=fread(fid,size,'float32')
end