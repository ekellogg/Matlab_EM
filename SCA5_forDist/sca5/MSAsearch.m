function [strseqnum,ats,best_align,alignment_trunc]=MSAsearch(pdb,chainID_trunc,alignment,truncate_or_not);
% usage: [seqnum,ats,best_align]=MSAsearch(pdb, 'A', aln); 
%        [seqnum,ats,best_align,aln_trunc]=MSAsearch(pdb, 'A',aln,1);
% 
% This function makes pairwise alignments between a query sequence (from
% the pdb file and chain ID) and every sequence in an MSA (alignment),
% finds the tophit sequence, and then attempts to make a residue number
% list (ats) that relates alignment numbering to structure numbering.  The
% function expects certain fields to exist in the pdb file that comprise
% the standard pdb format according to the pdb.org.  The function is therefore
% NOT guaranteed to work for user-specific or other non-standard pdb
% formats.
% 
% The default usage of this function is to read in a pdb structure file
% using getpdb.m or readpdb.m (functions in the biofinformatics toolbox)
% and then just make the residue number list (ats).  An alternative usage
% (by passing truncate_or_not=1) is to both make the residue number list
% (ats) and to truncate the MSA to the pdb sequence (returned in
% alignment_trunc).
%
% Arguments:
%
% Takes in:
%           (1) pdb, in pdb.org format (2) chainID, the particular sequence
%           to be used for truncation (3) alignment, the MSA (4)
%           truncate_or_not...1 to truncate to pdb, 0 to not.
% Returns:
%           (1) strseqnum, the MSA sequence number used for truncation (2)
%           ats, the residue position numbering according to the query
%           sequence, (3) best_align, the tophit pairwise alignment (4)
%           alignment_trunc, the truncated alignment if relevant.
%
%
% Authors: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%          Rohit Sharma
%          Kimberly Reynolds (kimberly.reynolds@utsouthwestern.edu)
%
% modified 1/2010 by R.R
% modified 2/2011 by Kim Reynolds
%
% Copyright R.Ranganathan 1999-2011
%**************************************************************************

% preliminaries
[x,y]=size(alignment);
if nargin <4
    truncate_or_not=0;
end

AA(1,:)={'ALA' 'ARG'  'ASN'  'ASP'  'CYS'  'GLN'  'GLU'  'GLY'  'HIS'  'ILE'...
    'LEU'  'LYS'  'MET'  'PHE'  'PRO'  'SER'  'THR'  'TRP'  'TYR'  'VAL'  'ASX'...
    'GLX'  'XAA'  'END' 'GAP'};
AA(2,:)={'A' 'R'  'N'  'D'  'C'  'Q'  'E'  'G'  'H'  'I'  'L'  'K'  'M'  'F'...
    'P'  'S'  'T'  'W'  'Y'  'V'  'B'  'Z'  'X'  '*' '-'};

%extracts relevant stuff from pdb. Makes resnumlist as cell array of
%strings and makes the amino acid sequence from the structure to be used
%for alignment querying.  The try-catch block tries to prevent one sort of
%exception...the case in which the pdb file has the filed pdb.Atom rather
%than the more usual pdb.Model.Atom.  The program takes into account the
%possibility that residue numbering includes special extra characters such
%as letters (e.g. residue 221A) by using the "iCode" field.
try
    if isfield(pdb.Model.Atom,'iCode')
        atom_fields=squeeze(struct2cell(pdb.Model.Atom));
        chainID=char(atom_fields(5,:));
        atomnum_fields=atom_fields(6:7,chainID==chainID_trunc);
        for i=1:size(atomnum_fields,2)
            atomnum_cat{i,:}=[num2str(cell2mat(atomnum_fields(1,i))) cell2mat(atomnum_fields(2,i))];
        end
        [resnumlist,indsort]=unique(atomnum_cat,'first');
        [aas]=atom_fields(4,indsort);
        [b,ix]=sort(indsort);
        resnumlist=resnumlist(ix);
        aa_model=aas(ix);
    else
        atom_fields=squeeze(struct2cell(pdb.Model.Atom));
        chainID=char(atom_fields(5,:));
        atomnum_field=atom_fields(6,chainID==chainID_trunc);
        [resnumlist,indsort]=unique(atomnum_field,'first');
        [aas]=atom_fields(4,indsort);
        [b,ix]=sort(indsort);
        resnumlist=resnumlist(ix);
        aa_model=aas(ix);
    end
catch
     if isfield(pdb.Atom,'iCode')
        atom_fields=squeeze(struct2cell(pdb.Atom));
        chainID=char(atom_fields(5,:));
        atomnum_fields=atom_fields(6:7,chainID==chainID_trunc);
        for i=1:size(atomnum_fields,2)
            atomnum_cat{i,:}=[num2str(cell2mat(atomnum_fields(1,i))) cell2mat(atomnum_fields(2,i))];
        end
        [resnumlist,indsort]=unique(atomnum_cat,'first');
        [aas]=atom_fields(4,indsort);
        [b,ix]=sort(indsort);
        resnumlist=resnumlist(ix);
        aa_model=aas(ix);
    else
        atom_fields=squeeze(struct2cell(pdb.Atom));
        chainID=char(atom_fields(5,:));
        atomnum_field=atom_fields(6,chainID==chainID_trunc);
        [resnumlist,indsort]=unique(atomnum_field,'first');
        [aas]=atom_fields(4,indsort);
        [b,ix]=sort(indsort);
        resnumlist=resnumlist(ix);
        aa_model=aas(ix);
     end
end

for i=1:numel(aa_model);aam_tmp(i)=AA(2,strmatch(aa_model(i),AA(1,:),'exact'));end
aa_model=cell2mat(aam_tmp);


% score each alignment sequence (removed of gaps) to the pdb sequence and
% find the tophit sequence to make the residue number list.  The typical
% expectation is that the pdb sequence is actually in the aligment and so
% we will find a 100% (or nearly 100%) match.
scores = zeros(1,size(alignment,1)); 
for rownum = 1:size(alignment,1)
    [scores(rownum),junk] = swalign(alignment(rownum,find(isletter(alignment(rownum,:)))), aa_model(isletter(aa_model)));
end
strseqnum = find(scores == max(scores));
if size(strseqnum,2) > 1
    strseqnum = strseqnum(1);
end

% Now we make the ats either with or without alignment truncation
    
if truncate_or_not==1  
    alignment_trunc = alignment(:,find(isletter(alignment(strseqnum,:))));
    disp(['     truncated alignment using sequence #' num2str(strseqnum) '  (score: ' num2str(max(scores)) ')']);
    [topscore, best_align, startat] = swalign(alignment_trunc(strseqnum,:), aa_model(isletter(aa_model)));
    best_align
    disp([' size of best_align is:  ' num2str(size(best_align,2))]);
    k=0;
    ats=cell(1,size(best_align,2));
    ats(1)=resnumlist(startat(2));
    for i=2:size(best_align,2);
        if best_align(2,i)=='|'|best_align(2,i)==':'
            k=k+1;
            ats(i)=resnumlist(startat(2)+k);
        elseif best_align(2,i)==' '& isletter(best_align(3,i))
            k=k+1;
            ats(i)=resnumlist(startat(2)+k);
        end
    end
        
else                 
    [topscore, best_align, startat] = swalign(alignment(strseqnum,:), aa_model(isletter(aa_model)));
    disp(['tophit is sequence #' num2str(strseqnum) '  (score: ' num2str(max(scores)) ')']);
    best_align
    disp([' size of best_align is:  ' num2str(size(best_align,2))]);
    k=0;
    ats=cell(1,size(best_align,2));
    ats(1)=resnumlist(startat(2));
    for i=2:size(best_align,2);
        if best_align(2,i)=='|'|best_align(2,i)==':'
            k=k+1;
            ats(i)=resnumlist(startat(2)+k);
        elseif best_align(2,i)==' '& isletter(best_align(3,i))
            k=k+1;
            ats(i)=resnumlist(startat(2)+k);
        end
    end
        
end

% now fix inconsistencies between the MSA and pairwise alignment that are
% inevitable since the ats is made from a pairwise alignment of the
% selected sequence from the MSA and the pdb sequence.  We do this by aia
% slightly bizarre way.  We make a profile of the pairwise alignemnt of the
% selected sequence from the MSA and the pdb model, and make a profile of
% the selected sequence as it exists in the MSA, and then make a profile to
% profile alignment.  
pf_1=seqprofile([best_align(1,:);best_align(3,:)]);
pf_2=seqprofile([alignment(strseqnum,:)]);
[prof,h1,h2]=profalign(pf_1,pf_2);
%old version:
%ats = ats(h2);

%ok... now - lets make a new ats.
%positions present in h2 (the alignment sequence), but not h1,
%should get assigned as gaps in the ats.
atsnew = cell(1,size(alignment,2));
atsnew(find(~ismember(h2,h1))) = cellstr(' ');
%things which are present in both h1 and h2 should get assigned from the
%previously constructed ats.
[isec, ia, ib] = intersect(h1,h2);
atsnew(ib) = ats(ia);
%note that things in h1 (the pdb sequence), but not in h2, should be left
%out of the ats, since they are (presumably) positions not
%present/truncated out of the alignment.
ats=atsnew;
    