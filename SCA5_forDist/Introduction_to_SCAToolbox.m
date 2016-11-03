%% SCA 5.0 MATLAB Toolbox

% This toolbox comprises a set of MATLAB scripts that include the main
% functions for executing the statistical coupling analysis (or SCA), a
% method for analyzing the pattern of evolutionary constraints in protein
% multiple sequence alignments (MSAs).  

% Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)
% Olivier Rivoire  (olivier.rivoire@ujf-grenoble.fr)

% August 2011

% Toolbox dependencies outside of basic MATLAB used in the tutorials:

% * Image Processing Matlab toolbox - for displaying correlation matrices. 
% 
% * Bioinformatics Matlab toolbox - for reading PDB files, and for using
% MSAsearch.m. 
%
% * Statistics Matlab toolbox - for some curve fitting operations used in
% the tutorials.
%
% * PyMOL molecular viewer - for displaying 3d structures.

%% An important note:  PLEASE READ!!

% Protein structure and function depend on functional interactions between
% amino acids.  The SCA approach rests on the idea that the relevant
% interactions might be revealed by correlations in the conservation of
% pairs of sequence positions in a representative sampling of members of a
% protein family (a multiple sequence alignment, or MSA). To perform this
% analysis, this toolbox comprises functions for: (1) calculating the
% degree of conservation of each position in a MSA taken independently (the
% first-order conservation statistic); (2) calculating a matrix of
% correlations in the conservation of pairs of positions (the second-order
% conservation statistic); (3) analyzing and interpreting the content of
% this matrix of both positional and sequence correlations in the MSA to
% deduce the statistically significant and functional patterns of amino
% acid correlations (termed "sectors", Halabi et al. (2009), Cell
% 138,774-786).  Work to date suggests that these sectors are novel
% features of protein three-dimensional structures, comprising physically
% contiguous networks of amino acids within proteins that contribute to
% core aspects of protein stability and function.

% SCA version 5 introduces some new general methods for defining and
% interpreting sectors in proteins.  The new work involves methods for
% mapping the relationship between patterns of correlation between amino
% acid positions (that define sectors) and patterns of correlation between
% sequence (that define functional or phylogenetic subfamilies). To explain
% the usage of SCA, we include four case studies that we feel are likely to
% represent the typical diversity of problems in sector identification in
% most proteins. However, we note that these applications are actively
% under study and development, and it would be inappropriate to treat SCA
% for now as a "black box" for sequence analysis.  We strongly suggest that
% users interested in applying these methods for new proteins first
% carefully go through the arguments and examples contained in the
% accompanying tutorials and technical note. Further information and
% discussions will be available on our laboratory web site
% (http://www.hhmi.swmed.edu/Labs/rr) as our studies develop.

%% Tutorials

% In this toolbox, we provide tutorials that describe application of SCA
% and sector identification in four different protein families (PDZ, the
% S1A serine proteases,  the Hsp70/110 family of molecular chaperones, and
% the G protein family of allosteric switches). Each case illustrates a
% distinct case of sector identification.  The tutorials are:

%   (1) The PDZ domain family of protein interaction modules
%       ("Tutorial_PDZ.m").  This tutorial illustrates the case of a
%       protein family with a single sector which appears to be a global
%       feature of members of the family.
%   (2) The S1A family of serine proteaes (trypsin, chmyotrypsin, etc.).
%       This tutorial illustrates the case in which multiple
%       quasi-independent sectors exist, each of which appears to control a
%       distinct functional aspect of these enyzmes. Like in PDZ, the
%       sectors seem to not be associated with obvious subfamilies of
%       proteins within the alignment.
%   (3) The Hsp70/110 family of molecular chaperones.  This family is known
%       to contain at least two distinct classes of sequences - the
%       allosteric Hsp70's and the non-allosteric Hsp110's - which
%       introduces the first instance of a highly inhomogeneous and
%       structured sequence space. This tutorial shows how we can take
%       advantage of this functionally relevant sequence divergence to find
%       the group of co-evolving positions (a sector) that corresponds to
%       the mode of sequence correlations that separates the Hsp70 and
%       Hsp110 subfamilies.  This sector seems to underlie the allosteric
%       mechanisn in the Hsp70 proteins.
%   (4) The G protein family of allosteric switches.  Is every instance of
%       sequence inhomogeneity associated with sectors?  In this tutorial,
%       we examine the case of the G protein superfamily, which contains
%       several clearly distinctly separated subfamilies...the Ras/Rho/Rab
%       family, the heterotrimeric G proteins, the translational elongation
%       and initiation factors, and the ADP-ribosylation factors.  Despite
%       the sequence inhomogeneity, we find a single sector in this protein
%       family that seems to be associated with the core shared function of
%       the G proteins - nucleotide dependent allosteric switching.


%% Notes

% A methodological supplement that reviews the basic formulae for SCA
% calculations is provided in the document entitled "Note 109 : Summary of
% SCA calculations". In addition to this note, information about the
% definition and identification of protein sectors can be found in the text
% and supplemental data of "Halabi et al., Cell 138, 774-786 (2009)" as
% well as "Smock et al., Mol. Syst. Biol. 6, 414 (2010)".

%% Inputs - Folder Contents

% This folder contains the alignments and structures that serve as inputs
% for the scripts. 

% Alignments in .free format:

% al_pdz.free      : PDZ alignment (240 sequences).
% al_S1A_1388.free : Serine protease alignment (1388 sequences).
% annot_S1A_1388   : Matlab file containing an annotation of the serine
%                    protease alignment
% al_hsp70.fasta   : Hsp70/110 alignment (926 sequences)
% ats_hsp70        : Position labels numbered according to E.coli DnaK, to
%                    be consistent with the associated structural model.
% G_alignment.free : The G protein alignment (678 sequences)
% G_labels.mat     : A MATLAB file with the names of sequences in the
%                    alignment.

% 3D structures from the PDB database:

% 1BE9.pdb           : The structure of PSD95-pdz3
% 3TGI.pdb           : rat trypsin in complex with BPTI
% DnaK-Sse1Model.peb : A model for the ATP-bound state of the E.coli DnaK
%                      Hsp70.  See Smock et al. for further details.
% 1Q21.pdb           : The structure of H-Ras in complex with GDP
% 5P21.pdb           : The structure of H-Ras in complex with GppNp

%% Functions - in the "sca5" folder

% See file headers for further descriptions of the following functions, 
% including the presentation of input/output options not used in the 
% Tutorials.

% Preprocessing tools (manipulation of alignments):
%
% get_seqs        : imports alignment in .free format.
% MSAsearch       : matches a query sequence to the alignment, with the option of
%                   truncating the alignment to the query sequence.
% lett2num        : converts the alignment in a numerical array.
% alnid           : truncates alignments to a user-specified max percent
%                   identity.

% Computational tools (calculations of conservation and correlations):
%
% sim_seq         : computes matrices of sequence correlations.
% cons            : computes positional conservation.
% sca5            : computes positional and sequence correlations ("SCA matrices").
% spectral_decomp : computes the spectral decomposition of the SCA
%                   positional correlation matrix.
% eigenvect       : computes the eigenvectors of a matrix.
% basic_ica       : an implementation of the independent component analysis

% Post-processing tools (graphical representation of the results):
%
% show_vect       : 1d representation of sectors based on 1 vector.  
% make_fig        : 2d representation of sectors based on 2 vectors.
% drawSS          : plots the organization of sectors by primary and
%                   secondary structure.
% write_pmlscript : writes a PyMol script for displaying the sectors in
%                   the three-dimensional structure.
% SCAcluster      : hierarchical clustering of the SCA correlation matrix.
% SCAhist         : histogram with fit for the SCA correlation matrix.
% scacursor       : graphical tool for examining SCA correlation matrix

% Accessory computational tools (for internal or alternative usages):
%
% binrep          : converts alignments to 3D binary tensors (for internal usage).
% weight_aln      : creates the weighted 3D alignment tensor. (internal
%                   usage)
% project_aln     : creates the 2D projected weighted alignment matrix
%                  (internal usage)
% DerivEntropy    : default weighting function (for internal usage).
% Entropy         : defines a weighting function based on Renyi entropies
%                   (for implementation of variations of the methods).

%% Outputs - Folder Contents

% Since they can be large and too cumbersome for distribution, this folder
% does not contain the files generated by executing the Tutorials.  Of
% course, running the tutorials will generate these files.

% Matlab workspaces associated with each tutorial:
% 
% work_pdz.mat   : outputs of 'Tutorial_pdz'
% work_sprot.mat : outputs of 'Tutorial_sprot'
% work_hsp70.mat : outputs of 'Tutorial_hsp70'
% work_G.mat     : outputs of 'Tutorial_G'

% PyMol scripts:
%
% sectors_pdz.pml   : The sector of the PDZ family on 1BE9.pdb
% sectors_sprot.pml : 3 sectors of the Serine Protease family on 3TGI.pdb

% sectors in hsp70 and G are not presented as pymol scripts.

