function [X3d] = binrep(algn)
% usage: [X3d]=binrep(algn); 
%
% This function converts a standard ASCII multiple sequuence alignment (M
% sequences x L positions) to a M x L x 20 3D binary tensor (X3d)
% representation.  x3d(n,i,a) = 1 if sequence n has amino acid a at
% position i, and is equal to 0 otherwise.
%
% This function is called in sca5.m.  
%
% ***********************
% Authors: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%          Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 8/2011
%
% Copyright R.Ranganathan 1999-2011
%**************************************************************************

% preliminaries
[N_seq,N_pos]=size(algn);
N_aa=20;

% lett2num.m converts an ascii alignment to numeric.  See file header of
% that program for details.

alg=lett2num(algn);
X3d=zeros(N_seq,N_pos,N_aa);
for a=1:N_aa
    X3d(:,:,a)=(alg==a);
end

end

