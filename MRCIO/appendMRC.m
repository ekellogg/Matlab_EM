function status= appendMRC(a, apix, fname)
% type is float or short
[fid,message]=fopen(fname,'r');
if fid == -1
    error('can''t open file');
    a= -1;
    return;
end
%read and replace header information
nx = fread(fid,1,'long');
ny = fread(fid,1,'long');
nz = fread(fid,1,'long');

if( nx ~= size(a,1) || ny ~= size(a,2) )
    error('sizes of files don''t match!')
    return;
end
type = fread(fid,1,'long');
fclose(fid);

fid = WriteMRCHeader(a,apix,fname,nz+size(a,3));
fclose(fid);

[fid,message] = fopen(fname,'r+');

sz = 0;
if(type == 1)
    sz = 2;
end
if(type == 2)
    sz = 4;
end
fseek(fid,0,'eof');
status=fwrite(fid,a,'float32');
fclose(fid);
end