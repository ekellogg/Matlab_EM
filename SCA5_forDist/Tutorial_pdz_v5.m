% SCA 5.0 Tutorial 1- Analysis of an alignment of the PDZ domain family.

% Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)
% Olivier Rivoire  (olivier.rivoire@ujf-grenoble.fr)

% - Aug. 2011
 
% This tutorial provides a practical walkthrough of the usage of the SCA
% Toolbox 5.0 using a multiple sequence alignment of the PDZ domain family
% as an example. The tutorial is in "cell" mode in MATLAB in which the
% commands corresponding to each section (or cell) can be executed by
% selecting the cell and clicking on the MATLAB command for "evaluate
% current cell". See Note 109 for the mathematical principles and details of
% the method.

% In general, the goal in SCA is to quantitatively deduce the information
% content of protein sequences.  The method is based on two ideas: (1) that
% the structural and functional properties fundamentally arise from the
% cooperative action of amino acids and (2) that the relative importance of
% interactions is indicated by their invariance (or conservation) in a
% representative statistical sampling of homologous sequences. Accordingly,
% this method attempts to deduce collectively evolving groups of conserved
% amino acid positions, termed "sectors".

% An important aspect of sector interpretation is to examine how sectors
% relate to subfamilies within the overall protein family; indeed, is a
% sector a global propety of the family or a property of functional
% specializations within subfamilies?  To address this, SCA v5.0 presents a
% general analysis that provides a quantitative mapping between
% correlations between positions (that define sectors) and correlations
% between sequences (that define subfamilies).

% This tutorial is the first in a series intented to illustrate different
% case studies of SCA calculation and intepretation.  This document will
% introduce many basic functions for completeness that will be assumed in
% subsequent tutorials; thus, we highly recommened going through this
% tutorial before the others.

addpath sca5
clear; close all

%% Step 1. Alignment loading and conditioning

% The toolbox contains a sample multiple sequence alignment of the PDZ
% domain family (240 sequences) in .free format. Load it using 'get_seqs'.
% Other file formats (e.g. FASTA or MSF) can be loaded using standard
% functions in the Bioinformatics Toolbox.

[labels_seq,algn_full]=get_seqs('Inputs/al_pdz.free');
N_seq=size(algn_full,1);

% A practical strategy is to truncate alignments to sequence positions with
% gap frequency no greater than 20%. This prevents trivial
% over-representation of gaps in the alignment and ensures that
% calculations are only made at largely non-gapped sequences positions.
% Other approaches to alignment truncation are possible that reflect this
% same general approach (see e.g. Tutorial_sprot).

cut_off=.2;
fraq_gaps=sum(isletter(algn_full)==0)/N_seq;
algn=algn_full(:,fraq_gaps<cut_off);
N_pos=size(algn,2);

% [OPTIONAL] - As part of the conditioning of the alignment, a good
% practice might be to eliminate redundant sequences prior to the analysis,
% so as to have no pair of sequences with a fraction of amino acids in
% common greater than some threshold (we often use a cutoff 90%).  However,
% the analysis seems quite robust to a small fraction of redundant
% sequences; for example, in this example we ignore this additional
% criterion).

% Labelling of positions / alternative method of truncation.
% It is often of value to have sequence positions numbered according to a
% specific member of the protein family (rather than per alignment
% numbering) to facilitate mapping the correlation data to the protein
% structure.  Given a structure file, first load the pdb file using the
% getpdb.m or pdbread.m functions, and make the position list using the
% function MSAsearch.m.  Note that MSAsearch expects a standard pdb file in
% the format specified by pdb.org.  This program is NOT guaranteed to work
% with user-specific or otherwise non-standard pdb files.

pdb_id='1BE9'; chain='A';
pdb=pdbread(['Inputs/' pdb_id '.pdb']);
[strseqnum,ats,best_align]=MSAsearch(pdb,chain,algn);

% 'strseqnum' identifies the alignment sequence that most closely matches to
% the query molecule (here chainID 'A' of pdb_1be9), 'ats' represents the
% alignment positions in pdb_1be9 numbering, and 'best_align' shows the
% pairwise alignment between the pdb sequence and the top-hit in the
% alignment. Note that MSAsearch can also be used to truncate the
% alignment to positions contained in the query molecule (see file header
% for usage). This represents another potential strategy for alignment
% truncation.

%% Step 2. Sequence similarity matrix

% The extent to which the alignment represents a "uniform" or "homogeneous"
% sampling of the sequence space can be judged from the structure of the
% correlations (or distances) between sequences. In one simple limit, the
% sequences would be equally dissimilar with no distinct cluster(s) of
% correlated sequences. In this case, sector identification amounts to
% examining the pattern of correlations between positions in the top
% eigenmodes of the SCA positional correlation matrix (as in Halabi et
% al.(Cell (2009), 138, 774-786), and see below).  This is the case for the
% PDZ domain family in this tutorial, the S1A serine protease family
% ("Tutorial_S1A.m"), and the DHFR family of metabolic enzymes
% ("Tutorial_DHFR.m").

% If the sequence space defined by the alignment shows strongly clustered
% patterns of sequences (that is, "inhomogeneous") then sector
% identification and interpretation can also be guided by the patterns of
% sequence divergence (as in Smock et al. MSB (2010), 6: 414). % In
% general, inhomogeneities in sequence divergence can arise due to either
% functionally or purely phylogenetic relationships between sequences in
% subfamilies.  A comparison between sector mappings and sequence
% correlations can in principle help distinguish the two.  Examples of such
% cases are given in the Hsp70 ("Tutorial_hsp70.m") and G protein
% ("Tutorial_G.m") families.

% To begin with a naive examination of the sequence space, we use the
% function sim_seq.m to compute a matrix (S) of similarity between pairs of
% sequences, such that S(s1,s2) gives the fraction of amino acids that are
% common between the sequences s1 and s2. Alternatively, we can define a
% matrix (G) of covariance between sequences (as in Halabi et al. Cell
% (2009), 138, 774-786), which is also available as optional output from
% simseq.m. These two matrices are vritually the same, differing only by a
% scaling factor.

[S]=sim_seq(algn);

% We make a histogram of the similarities between pairs of sequences (only
% half of the matrix S needs to be evaluated and the diagonal should be
% ignored):

listS=nonzeros(triu(S,1));
h_seqsim=figure; clf; 
set(h_seqsim,'Units','normalized','Position',[0 0.3 0.9 0.5],'Name','Sequence Correlations: PDZ');
subplot(1,2,1);hist(listS,N_pos/2);
xlabel('Pairwise SeqID','FontSize',14,'FontWeight','bold'); 
ylabel('number','FontSize',14,'FontWeight','bold'); grid on

% The histogram shows a reasonably narrow distribution with a mean pairwise
% identity between sequences of about 22% and a range of 10 to 40%,
% suggesting that most sequences are about equally dissimilar from other.
% We can further examine this assertion by direct visualization of the
% sequence similarity (S) matrix:

figure(h_seqsim); 
subplot(1,2,2); imshow(S,[0 1],'InitialMagnification','fit'); colormap(jet); colorbar;
title('SeqID', 'FontSize',12,'FontWeight','bold');

% The similarly matrix basically recapitulates what the histogram tell us;
% there are a few small clades of more related PDZ sequences, but in
% general, the alignment is comprised of a diverse and largely
% homogeneously diverged ensemble of sequences.  We will make this point
% again more meaningfully below.

%% Step 3. Positional conservation.

% In the current implementation, the degree of conservation of the
% different positions is measured in SCA by an information theoretic
% quantity called the Kullback-Leibler relative entropy - D(f(a,i)||q(a))
% See Note 103 for further details. In essence, this captures the
% divergence of the observed frequencies of amino acids at each position
% (f(a,i)) from their background frequencies in the non-redundant database
% of proteins (q(a)). The function cons.m makes this calculation, returning
% a 20 X N_pos matrix of relative entropies of each amino acid at each
% position in the multiple sequence alignment.  A global measure of
% conservation is similarly measured by D(f(i)||q') (where q' takes into
% account gaps, see Note 103, Section III.D), and the function cons.m
% returns this calculation in a 1 X N_pos vector of positional
% conservation.  cons.m also includes options to use a more general
% definition of entropy for conservation (called the Renyi entropy) and in
% principle to permit the use of user-defined functions. The file header of
% cons.m contains the relevant details.

% For the purpose of this tutorial, we only return the global
% Kullback-Leibler relative entropy of all amino acids at each position
% (Dglo, a vector of 1 X N_pos).

[D_glo]=cons(algn);

h_D=figure; set(h_D,'Units','normalized','Position',[0 0.6 0.4 0.4],'Name','Positional Conservation');clf
subplot(2,1,1);hist(D_glo,25); grid on;
xlabel('D (conservation)','FontSize',10,'FontWeight','bold'); 
ylabel('number','FontSize',10,'FontWeight','bold');
subplot(2,1,2);bar([1:numel(ats)],D_glo,'k'); grid on;
axis([0 numel(ats)+1 0 4]);
set(gca,'XTick',[1:10:numel(ats)]);
set(gca,'XTickLabel',ats([1:10:numel(ats)]));
xlabel('position (1BE9 numbering)','FontSize',10,'FontWeight','bold');
ylabel('D_i (conservation)','FontSize',10,'FontWeight','bold');

% The plots show a histogram of positional conservation values, and a bar
% graph of conservation values for each position.

%% Step 4. SCA calculations

% SCA v5.0 takes a multiple sequence alignment of a protein family as input
% and returns two main outputs: (1) a positional correlation matrix (Cp),
% which quantitatively indicates the correlated evolution of all pairs of
% positions in the alignment, and (2) a sequence correlation matrix (Cs),
% which indicates the pattern of similarity between all pairs of sequences.
% In SCA, both of these correlation matrices are weighted by the degree of
% conservation of amino acid positions.  Thus, the Cp matrix contains
% information about conserved correlations between pairs of positions, and
% the Cs matrix contains information about the similarity between pairs of
% PDZ sequences biased towards more conserved positions. See Note 109 for
% the mathematical principles and practical implementation of the method.

% The function sca5.m computes all this, and returns a MATLAB structure
% which contains all the relevant variables - the original alignment, the
% SCA positional correlation matrix Cp, the SCA sequence correaltion matrix
% Cs, and various important intermediate variables generated in the course
% of SCA calculation (the binarized, weighted, and projected alignments,
% the projection  matrix, the weight matrix, the weighting function, and
% optional parameter if any).  These variables are described in the
% technical tutorial to SCA calculations, and in the file header for sca5.m

[pdzsca]=sca5(algn);

%% Step 5. Spectral (or eigenvalue) decomposition

% To analyze the SCA positional correlation matrix, we carry out a
% mathematical technique called spectral (or eigenvalue) decomposition. The
% motivation is to relaize that the existence of non-trivial correlations
% between positions indicates that treating the amino acids as the basic
% units of proteins is not the most informative representation. Instead, we
% should seek a reparameterization of the protein in which the units of
% proteins are the collective groups of amino acids that coevolve per the
% positional correlation matrix. These collective groups of coevolving
% residues are called "sectors", and are proposed to be the fundamental
% units in the evolutionary "design" of proteins (Halabi et al., Cell
% (2009), 138, 774-786).

% Eigenvalue decomposition is the simplest first step in achieving this
% reparamterization (see "Tutorial_S1A.m" for a more sophisticated
% approach).  This decomposition is always possible for any square positive
% semi-definite matrix (such as a correlation matrix) and mathematically
% transforms a current representation of a system in which variables are
% correlated (for example, the sequence positions) into new variables that
% now have the property of being uncorrelated to each other. These new
% variables are linear combinations of the original variables (in our case,
% groups of sequence positions) and represent the more informative
% parameteriation of the system. In addition the eigenvalue decomposition
% provides a basis to order the new transformed variables by magnitude of
% capturing the information content of the original corrleation matrix.

% How deos eigenvalue decomposition work? The original matrix is written as
% a product of three matrices: X=VDV', where D is a diagonal matrix of
% so-called eigenvalues and columns of V contain the associated so-called
% eigenvectors.  The eigenvectors contain the weights for linearly
% combining the original variables into each new transformed variable and
% the eigenvalue associated with each eigenvector indicates the magnitude
% of information in the original matrix captured. In other words, for the
% SCA positional correlation matrix, each eigenvector represents a weighted
% combination of sequence positions (an "eigenmode"), and the
% associated eigenvalue indicates the statistical importance of that mode.

% How many eigenvalues are statistically significant?  For alignments in
% which the number of sequences is not large compared to the number of
% positions (the usual case), we expect that most eigenvalues are simply
% explained by statistical noise coming from finite sampling of sequences.
% To check this, we compare the spectral decomposition for the actual
% alignment with that for many instances of randomized alignments in which
% amino acids are scrambled independently down each column.  This
% manipulation removes all functional correlations and retains only the
% spurious correlations that are possible due to finite sampling.  The
% function spectral_decomp.m carries out this calcualtion and returns a
% structure with the eigenvalue decomposition of actual and randomized
% alignments. It also makes a plot of the eigenspectra.

[spect]=spectral_decomp(pdzsca,100);

% In the MATLAB structure spect, spect.lbdpos and spect.lbdrndpos contain the
% eigenvalues of the actual and randomized alignments, and spect.evpos and
% spect.evrndpos contains the corresponding eigenvectors. For the PDZ family,
% it is clear from the plotted histograms of eigenvalues that most
% eigenvalues are explained by just finite sampling limitations (the red
% trace is the random expecatation), and that just the top two or three
% eigenvalues are statistically significant.  These top eigenmodes can be
% examined for any patterns of positional correlations.

%% Step 6a.  Structure of top eigenmodes

% Sectors are empirically defined by examining the pattern of positional
% contributions to the top few eigenvectors. We will just examine the
% structure of the three top eigenmodes here.

% a 3-D plot of the top three eigenvectors
h_3Dtopmodes=figure; set(h_3Dtopmodes,'Units','normalized','Position',[0 0.7 0.3 0.4],'Name','Top Eigenmodes - 3D'); clf; 
scatter3(spect.evpos(:,1),spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,1)+.01,spect.evpos(i,2)+.01,spect.evpos(i,3)+.01,ats(i));end;hold off
az=136;el=20;view(az,el);
xlabel('ev 1','FontSize',12,'FontWeight','b');
ylabel('ev 2','FontSize',12,'FontWeight','b');
zlabel('ev 3','FontSize',12,'FontWeight','b');


% 2-D plots of the top three eigenvectors
h_2Dtopmodes=figure; set(h_2Dtopmodes,'Units','normalized','Position',[0 0 1.0 0.4],'Name','Top Eigenmodes-2D'); clf; 
subplot(1,3,1);
scatter(spect.evpos(:,1),spect.evpos(:,2),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,1)+.01,spect.evpos(i,2)+.01,ats(i));end;hold off;grid on
xlabel('ev 1','FontSize',12,'FontWeight','b');ylabel('ev 2','FontSize',12,'FontWeight','b');
subplot(1,3,2);
scatter(spect.evpos(:,1),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,1)+.01,spect.evpos(i,3)+.01,ats(i));end;hold off;grid on
xlabel('ev 1','FontSize',12,'FontWeight','b');ylabel('ev 3','FontSize',12,'FontWeight','b');
subplot(1,3,3);
scatter(spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,2)+.01,spect.evpos(i,3)+.01,ats(i));end;hold off;grid on
xlabel('ev 2','FontSize',12,'FontWeight','b');ylabel('ev 3','FontSize',12,'FontWeight','b');

% In the case of the PDZ domain, examination of the top modes shows three
% things: (1) most positions contribute weakly (weights near zero) and
% therefore are predicted to evolve nearly independently, (2) a few
% positions (~20%) contribute to the tails of the eigenvector
% distributions, and (3) the contributing positions form a scattered
% pattern of positional weights that do not obviously indicate patterns of
% multiple independent sectors.  Rather, the data are consistent with a
% single sector in which constituent positions are involved in a
% hierarchical pattern of correlations with each other.

%% Step 6b. Analysis by hierarchical clustering [OPTIONAL]

% We can also more intuitively see this by hierarchical clustering of the
% Cp matrix. In this approach, two positions will be juxtaposed if they
% show a similar profiles of correlation with other positions.  Clustering
% is carried out by SCAcluster.m, and a useful program -  scacursor.m - for
% interrogating clusters (by Bill Lane, Seth Darst's lab, The Rockefeller
% University)

[p_pos,l_pos,sort_pos,Csorted]=SCAcluster(pdzsca.Cp,ats,1);
figure(gcf);scacursor(ats(sort_pos),ats(sort_pos),Csorted);

% The figure number referenced prior to scacursor should contain the SCA
% correlation matrix, and it is necessary to pass the position labels as a
% cell array. Use the left and right arrow keys to adjust the color scale.
% As documented in its file header, scacursor.m contains capabilities for
% interactively extracting cluster positions in a manner suitable for cut
% and paste into the molecular graphics program PyMol (www.pymol.org).  If
% this function causes errors at runtime, comment it out in executing this
% cell.

% Again, clustering shows that the statistical structure of the Cp matrix
% seems to be dominated by one group of hierarchically correlated
% positions.

%% Step 6c. A mapping between positional and sequence correlations

% An aspect of SCA version 5 (sca5) is the ability to map between the space
% of position correlations (whose top eigenmodes define the sectors) and
% the space of sequence correlations (whose top eigenmodes tell us about
% the sub-family architecture, if any) in the alignment.  This mapping can
% be used generally to study how the amino acid composition of sectors is
% related to functional or phylogenetic differences between sequences.

% The key in making this mapping between sequence correlations and
% positional correlations is a generalization of eigenvalue decomposition
% known as singular value decomposition.  To explain, the general
% mathematical representation of a multiple sequence alignment is a
% 3D-binary tensor x3d (M sequences x L positions x 20 amino acids).
% x3d(n,i,a) = 1 if sequence n has amino acid a at position i, and is equal
% to 0 otherwise.   In SCA5, this 3D tensor is reduced (or "projected") to
% a 2D matrix of weighted values that represents each amino acid at each
% position. This project sequence alignment (pwX) is returned in the MATLAB
% structure output by sca5.m (here, pdzsca.pwX).

% By the SVD, the alignment matrix pwX (M sequences x L positions) can be
% written as a product of three matrices: pwX=U*sv*V', where U is a M x M
% matrix whose columns contain the eigenvectors of the sequence correlation
% matrix, and V is an L x L matrix whose columns contain the eigenvectors
% of the positional correlation matrix.  sv is a diagonal M x L matrix, and
% contains the so-called "singular values" that are related to the
% eigenvalues of the sequence and positional correlation matrices.  The SVD
% can be easily calculated in MATLAB:

[U,sv,V]=svd(pdzsca.pwX);

% The important concept is that the matrix Pi=U*V' provides a mathematical
% mapping between the positional and sequence correlation matrices: 

N_min=min(N_seq,N_pos);
Pi=U(:,1:N_min)*V(:,1:N_min)';
U_p=Pi*spect.evpos;

% Thus, if an eigenmode of the positional correlation matrix (a columns in
% spect.evpos) describes a sector (a group of coevolving amino acid
% positions), then the corresponding column of U_p will contain the the
% pattern of sequence divergence (if any) that underlies the sector
% definition.  We can make this mapping:

h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Seq Correlations by Positional Correlations'); clf; 

h_SectSeq(1)=subplot(1,2,1)
scatter3(spect.evpos(:,1),spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,1)+.01,spect.evpos(i,2)+.01,spect.evpos(i,3)+.01,ats(i));end;hold off
az=58;el=30;view(az,el);
xlabel('ev 1','FontSize',12,'FontWeight','b');
ylabel('ev 2','FontSize',12,'FontWeight','b');
zlabel('ev 3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2)
scatter3(U_p(:,1),U_p(:,2),U_p(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
az=58;el=30;view(az,el);
xlabel('Seq 1','FontSize',12,'FontWeight','b');
ylabel('Seq 2','FontSize',12,'FontWeight','b');
zlabel('Seq 3','FontSize',12,'FontWeight','b');

% For the PDZ domain family, this mapping is rather boring.  There is no
% particular structure to positional correlations (as mentioned above) and
% no structure in the sequence space of the protein family either. Instead,
% sequences just seem to homogeneously occupy the space spanned by the top
% few eigenmodes of the sequence correlation matrix. Thus, our conclusion
% is a single hierarchical sector that is a global property of this protein
% family.  More interesting cases of this sort of sector - sequence mapping
% are described in the other tutorials.

% Here, we mapped the top eigenmodes of sequence correlations onto the top
% eigenmodes of the positional correlations. But this mapping can be
% inverted; if an eigenmode of the sequence correlation matrix describes a
% separation of two functionally distinct classes (or subfamilies) of
% sequences, then the corresponding eigenmode of the positional correlation
% matrix will reveal the sector that is responsible for this sequence
% divergence.  We will use this mode of sector identiication in
% "Tutorial_hsp70.m", and see Smock et al. MSB 6,414.

%% Step 6d. Sector definition

% We will define the PDZ sector based on the distribution of positional
% weights to the top (first) eigenmode of the Cp matrix.  Empirical
% examination shows that this distribution is well-fit by the lognormal
% distribution. Note that no deep mechanistic basis for the use of this
% distribution is implied here; indeed, understanding the statistical
% process that generates the observed distribution of correlations in
% natural proteins is a major research goal.  Here, we use this fitting
% simply to provide a logical and systematic basis for defining sector
% positions.

% We fit top eigenvector to the lognormal distribution and define the
% sector by making a cumulative density function (cdf) from the fitted
% distribution and picking a cutoff in the tail (at 80%). So, this means we
% are choosing to represent the residues belonging to the top 20% of the
% cdf.  This cutoff illustrates the basic properties of sectors...a few
% positions in the tail of the distibution of eigenvector weights which (as
% we show further below) comprise distributed physically contiguous
% networks of atoms in the protein structure linking the main active site
% to distantly positioned surface sites. 

% figure and probability cutoff
h_secdef=figure;
set(h_secdef,'Units','normalized','Position',[0 1 .5 0.3],'Name','Top Eigenmode'); clf;
p_cutoff=0.8; % the cdf cutoff
secpos = [];

% histogram of the data
xhist=[0:.01:.4]; % make bins for the histogram based on the full range of the data
[yhist]=hist(spect.evpos(:,1),xhist); bar(xhist,yhist./N_pos,'k');hold on;grid on

% distribution fitting
pd=fitdist(spect.evpos(:,1),'lognormal');
x_dist=[min(xhist):(max(xhist)-min(xhist))/100:max(xhist)];
area_hist=N_pos*(xhist(2)-xhist(1)); % for proper scaling of the pdf
pdf_jnk=pdf(pd,x_dist);
scaled_pdf=area_hist.*pdf_jnk;
plot(x_dist,scaled_pdf./N_pos,'r-','LineWidth',1.5);

% here, we make the cdf, and define sectors:
cdf_jnk=cdf(pd,x_dist);
clear sec cutoff_ev
[jnk,x_dist_pos_right]=min(abs(cdf_jnk-(p_cutoff)));
cutoff_ev = x_dist(x_dist_pos_right)';
% we obtain the indices of sector positions given the cutoffs
[sec.def] = find(spect.evpos(:,1)>cutoff_ev);
sprintf('%g+',str2num(char(ats(sec.def))))
sec.cutoff=cutoff_ev;
figure(h_secdef); line([cutoff_ev cutoff_ev],[0 max(yhist)/N_pos],'LineWidth',1,'LineStyle','--','Color','b');
sec.col=2/3;
%sec.vec=spect.evpos(:,1)-cutoff_ev;

% the eigenvector distribution for randomized alignments, scaled for
% comparison by the ratio of eigenvalue magnitude for randomized and actual
% alignments.

clear tmp
N_samples=100;
for i=1:N_samples;tmp(i,:)=squeeze(spect.lbdposrnd(i,1).*spect.evposrnd(i,:,1));end
[yhist_rnd]=hist(tmp(:)./spect.lbdpos(1),xhist); 
plot(xhist,yhist_rnd/(N_samples*N_pos),'g','LineWidth',1.5);
hold off

%% Step 10. Structural analysis of sectors

% One of the central findings is that sectors are not obviously related to
% current classification of proteins by primary, secondary, and tertiary
% structure. Indeed, sectors comprise positions that can be widely
% distributed in primary structure, that do not follow known patterns of
% contacts due to secondary structure, and do not follow tertiary
% structural motifs and subdomain architectures. Instead, they appear to be
% new features of three-dimensional structure that follow only the rule
% that they are physically connected in the three-dimensional structure, at
% least in one physiological state of the protein. To see this in context
% of the PDZ domain, we map sector positions by primary, secondary, and
% tertiary structure representations.

% Organization by the primary/secondary structure:

h_PriSecStr=figure;
set(h_PriSecStr,'Name','Sectors: Primary and Secondary Structure','Position',[41 37 841 266]);
clf
subplot(2,1,1);
show_vect(D_glo,sec,ats);
pdzss.Sheet=[1 2 3 6 7]; pdzss.Helix=[2 3];
figure(h_PriSecStr); hsub=subplot(2,1,2);hold on;drawSS(pdb,ats,pdzss,sec);hold off

% Representation on a 3D structure: a script is generated that displays the
% sectors when executed with the software PyMol (www.pymol.org). Once PyMol
% is installed, you should just be able to double click on the
% "Outputs/sector_pdz.pml".

filename='Outputs/sector_pdz';
write_pmlscript(pdb_id,chain,ats,sec,filename);
% Here we display the ligand (chain B) in stick bonds and color it by
% element with carbons in yellow:
fid=fopen([filename '.pml'], 'a');
fprintf(fid,'set stick_radius, 0.4\n');
fprintf(fid,'show stick, (chain B)\n');
fprintf(fid,'util.cbay main and chain B\n');
% and, we remove the n-terminal and c-terminal extensions that are not part
% of the core PDZ domain
fprintf(fid,'remove main and chain A and resi 300-308\n');
fprintf(fid,'remove main and chain A and resi 394-end\n');
% some other cosmetics that we prefer
fprintf(fid,'set two_sided_lighting, 1\n');
fprintf(fid,'set sphere_quality, 2\n');
fprintf(fid,'set surface_quality, 2\n');
fprintf(fid,'set stick_quality, 10\n');
fprintf(fid,'set cartoon_oval_quality, 10\n');
fprintf(fid,'set cartoon_loop_quality, 6\n');
% the following sets the default view
fprintf(fid,'set_view (\\\n');
fprintf(fid,'-0.880436897,    0.134605691,    0.454657406,\\\n');
fprintf(fid,'0.111319803,    0.990739465,   -0.077750385,\\\n');
fprintf(fid,'-0.460912317,   -0.017842097,   -0.887267053,\\\n');
fprintf(fid,'0.000000000,    0.000000000, -104.157783508,\\\n');
fprintf(fid,'35.284568787,   61.732337952,   29.176891327,\\\n');
fprintf(fid,'88.666572571,  119.648986816,  -20.000000000 )\n');

fclose(fid);

% The program write_pmlscript.m is not particularly sophisticated or ideal
% in taking advantage of features in PyMol to effectively represent
% sectors.  A more powerful but manual way to examine sectors is simply to
% write out sector residues in a format suitable for cutting and pasting
% into PyMol and work interactively between MATLAB and PyMol.  For example:

sprintf('%g+',str2num(char(ats(sec.def))))

% This command produces a residue list that can be pasted into a selection
% or object definition in PyMol.  To gradient color sectors by strength of
% contribution (as in Fig.3 of Halabi et al. Cell (2009), 138, 774-786), we
% use the excellent web-downloadable Python scripts written by Dr. Robert
% Campbell of Queen's University, Canada
% (http://pldserver1.biochem.queensu.ca/~rlc/) to (1) write out the weights
% for sector residues into the B-factor column (data2Bfactor.py) and (2) to
% color by a user-specified gradient (color_b.py).

%% Step 11.  Conclusion

% In the PDZ domain family, we find a single hierarchically organized
% sector that comprises a physically contiguous network of atoms that links
% the ligand binding site to two distantly positioned surface sites.  There
% is no evidence that this sector emerges from any subfamily architecture
% in the alignment, and therefore we conclude that it is a global property
% of the PDZ family.

%% 

save Outputs/work_pdz