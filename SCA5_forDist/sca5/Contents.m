% The SCA 5.0 MATLAB Toolbox
% v5.0 Aug, 2011
%
% Authors:
%   Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)
%   Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
%
%
% This toolbox comprises a set of MATLAB scripts that include the main
% functions for executing the statistical coupling analysis (or SCA), a
% method for analyzing the pattern of evolutionary constraints in protein
% multiple sequence alignments (MSAs). 
%
% **********************************************************************
% NOTE:  It is strongly suggested that prior to using this toolbox, users
% read the material in the file "Introduction_to_SCAToolbox.m" and go
% through the tutorials that accompany this toolbox.
% **********************************************************************
%
% External dependencies (none necessary for the core functions): 
% * Image Processing Matlab toolbox - for displaying correlation matrices. 
% 
% * Bioinformatics Matlab toolbox - for reading PDB files, and for using
% MSAsearch.m. 
%
% * Statistics Matlab toolbox - for fitting distributions in analysis of
% the correlation matrix in the tutorials.
%
% * PyMOL molecular viewer - for displaying 3d structures.
%
%
% FUNCTIONS (click on each to see the file headers)
%
% Files
%
%   alnid           - usage: [aln90]=alnid(aln,0.9);;
%                            This function truncates alignments to a
%                            user-specified maximum identity betweeen
%                            sequences.
%   basic_ica       - usage  w=basic_ica(x,r,Niter);
%                            This function is a relatively simple
%                            implemenation of an algorithm for independent
%                            component analysis
%
%   binrep          - usage: [X3d]=binrep(algn);
%                            This function takes in an aligment in text
%                            format (algn) and returns the 3D binary tensor
%                            representation.
%
%   cons            - usage: [D,Dmat]=cons(algn);
%                            This function computes the conservation of
%                            amino acids at positions in a multiple
%                            sequence alignment (algn), and returns the
%                            data either as a vector (D) of global
%                            conservation of amino acids at each position,
%                            or a matrix (Dmat) of conservation values for
%                            each amino acid at each position.
%
%   DerivEntropy    - usage: D=DerivEntropy(f,q)
%                            Computes the derivative of relative entropy
%                            dD(f||q)/df; an interal function for sca_mat.
%
%   drawSS          - usage: drawSS(pdb,alignment_to_structure)
%                            This function takes in a pdb file and a list
%                            of positions (ats) and outputs a graph showing
%                            the secondary structure pattern of the protein
%                            trimmed to the list of selected positions.
%
%   eigenvect       - usage : [w,r]=eigenvect(A)
%                             Returns all (or the only the first k)
%                             principal eigenvectors of matrix A, w(:,1),
%                             w(:,2), ..., and their corresponding
%                             eigenvalues.
%
%   Entropy         - usages: D=Entropy(f,q)
%                             Relative entropy D(f|q) if a is not specified
%                             or a=1. Renyi (relative) entropy with
%                             parameter a otherwise.  The basic measure of
%                             conservation in SCA.
%
%   get_seqs        - usage: [seqID, alignment]=get_seqs(filename)
%                            Imports an alignment in .free format. In this
%                            format, an alignment is represented as
%                            follows: each line should contain a seqID, a
%                            tab character, the sequence comprised of the
%                            20 amino-acids and a gap denoted by a period
%                            or a dash.  Each line is separated by a
%                            paragraph mark. 
%
%   lett2num        - usage: msa_num=lett2num(msa_lett)
%                            Translates an alignment from a representation
%                            where the 20 natural amino acids are
%                            represented by letters to a representation
%                            where they are represented by the numbers
%                            1,...,20.
%
%   make_fig        - usage: []=make_fig_col(v1,v2,sec);
%                            2d scatter plot with optional gradients of
%                            color and size for the points. used to show
%                            patterns of positional correlation and
%                            sequence similarity.
%
%   MSAsearch       - usage: [seqnum,ats,best_align]=MSAsearch(pdb, 'A', algn); 
%                            used either for truncation of an alignment per
%                            residues corresponding to a specific atomic
%                            structure, or for construction of a position list
%                            relating alignment position to residues
%                            corresponding to a specific atomic structure.
%
%   project_aln     - usage: [pwX]=project_aln(algn,wX,W);
%                            This function takes in an alignment (algn),
%                            the weighted tensor representation of the
%                            alignment (wX, output by weight_aln), and the
%                            weight matrix (W, output by weight_aln), and
%                            produces the 2D weighted projected alignment
%                            matrix (pwX).
%
%   sca5            - usages: [SCA_struct, Csca]=sca5(algn);
%                             This function computes the SCA correlation
%                             matrices Cp and Cs and various intermediate
%                             variables given a protein multiple sequence
%                             alignment (algn). The output is returned in
%                             an MATLAB structure file.  sca5 calls binrep,
%                             weight_aln, and project_aln.  Optionally, it
%                             can also return the positional correlation
%                             matrix Csca, which corresponds to the matrix
%                             returned by earlier version of the SCA method
%                             (v3.0-4.5).
%           
%   SCAcluster      - usage: [p,l,sort_order,sorted]=SCAcluster(matrix,pos,1.0,jet,1);
%                            Two dimensional hierarchical clustering of SCA
%                            correlation matrix using city-block distance
%                            metric and compete linkage. 
%
%   scacursor       - usage: scacursor(pos_sorted,pos_sorted,C_sorted);
%                            a tool for examining the SCA correlation matrix, and
%                            if clustered, for extracting cluster composition in
%                            a manner suitable for cut and paste into PyMol, the
%                            molecular graphics program.
%
%   show_vect       - usage: []=show_vect(V,sec);
%                            representation of sectors by specified colors
%                            on the primary structure of the protien
%                            family.
%
%   sim_seq         - usage: [sim_mat,G]=sim_seq(algn)
%                            Returns matrices of correlation and similarity
%                            between sequences in a multiple sequence
%                            alignment.
%
%   spectral_decomp - usage: [spect]=spectral_decomp(SCA_struct,N_samples);
%                            This function takes in the SCA_struct
%                            (returned by sca5), and computes the spectrral
%                            (or eigenvalue) decomposition of the SCA
%                            positional correlation matrix and the number of
%                            trials of this decomposition for randomized
%                            alignments specified by N_samples.
%
%   weight_aln      - usage: [wX,W]=weight_aln(x3d);
%                            This function takes in the 3D binary tensor
%                            represenatation of the MSA (output by binrep)
%                            and returns the weighted 3D alignment tensor
%                            (wX) and the weight matrix (W).
%
%   write_pmlscript - usages: []=write_pmlscript(pdb_id,chain,labels_pos,sec,output)

%
