function [labels,seqs]=get_seqs(filename)
% usage: [seqID, alignment]=get_seqs(filename)
%
% Imports an alignment in .free format. In this format, an alignment is
% represented as follows: each line should contain a seqID, a tab character,
% the sequence comprised of the 20 amino-acids and a gap denoted by a
% period or a dash.  Each line is separated by a paragraph mark. For FASTA
% or other file formats, use the utilities provided in the MATLAB
% bioinformatics toolbox.
%
% Author: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
% 
%
%**************************************************************************

[labels,seqs]=(textread(filename,'%s%s'));
labels=char(labels);
seqs=char(seqs);