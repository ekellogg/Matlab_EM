function [new_aln,idkeeplist,idx]=alnid(aln,cutoff,idmat);
% usage: [aln90]=alnid(aln,0.9);
%        [aln90,ind_keep]=alnid(aln,0.9);
%**************************************************************************
%**************************************************************************
% Author: Bill Russ (russ@chop.swmed.edu)
%
% BR 02/13/2007, modified to operate from a pre-calculated identity matrix
% RR 2/2009, added comments
%
% This function takes in an alignment (aln), a fractional cutoff for
% identity (cutoff), and optionaly, a pre-calculated identity matrix
% (idmat) and outputs a truncated alignment (new_aln) in which every
% sequence is less than the cutoff top hit identity to any other sequence.
% Probably could/should be modified to minimize usage of loops by
% vectorization.
%
%
%**************************************************************************
%**************************************************************************
% preliminaries
if ischar(aln);
    nseqs = size(aln,1);
else
    if iscell(aln)'
        nseqs = numel(aln);
    end
end

idkeeplist = ones(1,nseqs);

% make list of sequences to keep
if nargin < 3
    for i = 1:nseqs-1
        for j = i+1:nseqs;
            s1 = aln(i,:);
            s2 = aln(j,:);
            h = (s1~='-') | (s2~='-'); % only avoid sites with double gaps
            f = 1-(sum( (s1~=s2) & h ) / sum(h));
            idx(i,j) = f;idx(j,i)=f;
            if (f > cutoff)
                idkeeplist(j)=0;
                j = i+1;
            end
        end
    end
else
    idx=idmat;
    for i = 1:nseqs-1
        for j = i+1:nseqs;
            f = idx(i,j);
            if (f > cutoff)
                idkeeplist(j)=0;
                j = i+1;
            end
        end
    end
end

% truncate alignment
if ischar(aln);
    new_aln = aln(find(idkeeplist),:);
else
    if iscell(aln);
        new_aln = aln(find(idkeeplist));
    end
end


sumidx = sum(idx);
xx = find(sumidx==0 | isnan(sumidx));
idkeeplist(xx)=0;
            
