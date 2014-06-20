function a= readMRCfile(fname,first_index,last_index)
[fid,message]=fopen(fname,'r');
if fid == -1
    error('cant open file');
    a= -1;
    return;
end
nx=fread(fid,1,'long');
ny=fread(fid,1,'long');
nz=fread(fid,1,'long');
type= fread(fid,1,'long');
fprintf(1,'nx= %d ny= %d nz= %d type= %d', nx, ny,nz,type);
if( nargin == 3 )
   startpos = (first_index*nx*ny);
   size = (last_index-first_index)*(nx*ny);
elseif( narg == 2 ) %read in contents starting from index first_index
   startpos = (first_index*nx*ny);
   size = nx*ny*nz;
else
   startpos = 0;
   size = nx*ny*nz;
end

% Shorts
if type== 1
    % Seek to start
    status=fseek(fid,1024+(startpos*2),-1); %seek to startpos in bytes
    a=fread(fid,size,'int16');
end
%floats
if type == 2
    % Seek to start
    status=fseek(fid,1024+(startpos*4),-1);
    a=fread(fid,size,'float32');
end
fclose( fid);
nimg = (size)/(nx*ny);
a= reshape(a, [nx ny nimg]);
if nimg == 1
    a= reshape(a, [nx ny]);
end
