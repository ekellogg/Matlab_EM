% SCA 5.0 Tutorial 4 - Analysis of an alignment of the small G protein
% family.

% Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)
% Olivier Rivoire  (olivier.rivoire@ujf-grenoble.fr)

% - Aug 2011

% This tutorial presents SCA for the G protein family.  This family
% illustrates an interesting case - sequences are distributed into distinct
% subfamilies, but this subfamily architecture is irrelevant to the pattern
% of positional correlations.  Indeed, there is one apparent sector that
% seems to be a global property of members of this family regardless of the
% subfamily architecture.  

% This tutorial builds on the concepts introduced in 'Tutorial_pdz',
% 'Tutorial_sprot' and 'Tutorial_Hsp70'.  As such, we suggest that these
% tutorials be examined carefully before working through this protein
% family.

addpath sca5
clear; close all

%% Step 1. Alignment loading 

% The toolbox contains a sample multiple sequence alignment of the G
% protein family (comprising 678 sequences) in .free format, which is
% already trimmed to have positions with greater than 20% gaps removed.
% Load it using 'get_seqs'.

[algn]=get_seqs('Inputs/G_alignment.free');
load Inputs/G_labels.mat
[N_seq,N_pos]=size(algn);

% We make a residue numbering list that maps alignment positions to the
% numbering system of the GDP state of ras (pdb 1Q21), the GTP state of
% cdc42 (1NF3), or the GTPgammaS bound state of Gialpha. These are
% alternative structures for examining the architecture of setors.  See the
% "Tutorial_pdz.m" for details on this step.

pdb_1q21=pdbread(['Inputs/1Q21.pdb']);
[strseqnum_1q21,ats_1q21,best_align_1q21]=MSAsearch(pdb_1q21,'A',algn);

%% Step 2. Sequence similarities

% To begin with a naive examination of the sequence space described by the
% alignment, we use the function sim_seq.m to compute a matrix (S) of
% similarity between pairs of sequences, such that S(s1,s2) gives the
% fraction of amino acids that are common between the sequences s1 and s2.
% Alternatively, we can define a matrix (G) of covariance between sequences
% (as in the Halabi et al. Cell (2009), 138, 774-786), which is also
% available as optional output from simseq.m. These two matrices are
% vritually the same, differing only by a scaling factor.

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

% The histogram shows a multi-peaked distribution with
% groups of sequence showing mean pairwise identities between sequences
% ranging from 10% to greater than 50%.  This pattern already suggests the
% existence of multiple subfamilies. We can further examine this assertion
% by direct visualization of the sequence similarity (S) matrix:

figure(h_seqsim); 
subplot(1,2,2); imshow(S,[0 1],'InitialMagnification','fit'); colormap(jet); colorbar;
title('SeqID', 'FontSize',12,'FontWeight','bold');

% The similarity matrix again indicates the presence of several subfamilies
% of sequences.  As shown in the case of the Hsp70/110 family, the rich
% subfamily structure in this alignment could in principle tell us more
% about specific sectors related to various aspects of functional
% divergence between subfamilies.  It is also possible that the subfamily
% architecture is purely phylogenetic in nature and unrelated to sectors.
% Indeed, here we show that despite the sequence inhomogenieties, the G
% protein family seems to contain one sector that is related to the core,
% shared function of members of this protein family - to act as
% nucleotide-dependent switches that control interaction with effector
% molecules at specific molecular surfaces.  

%% Step 3. Positional Conservation

% We compute the degree of conservation of each position by the
% Kullback-Leibler relative entropy - D(f(a,i)||q(a)) See Note 109 and
% Tutorial_pdz.m for further details.

[D_glo]=cons(algn);

h_D=figure; set(h_D,'Units','normalized','Position',[0 0.6 0.5 0.4],'Name','Positional Conservation');clf
subplot(2,1,1);hist(D_glo,25); grid on;
xlabel('D (conservation)','FontSize',10,'FontWeight','bold'); 
ylabel('number','FontSize',10,'FontWeight','bold');
subplot(2,1,2);bar([1:numel(ats_1q21)],D_glo,'k'); grid on;
axis([0 numel(ats_1q21)+1 0 4]);
set(gca,'XTick',[1:10:numel(ats_1q21)]);
set(gca,'XTickLabel',ats_1q21([1:10:numel(ats_1q21)]));
xlabel('position (1Q21 numbering)','FontSize',10,'FontWeight','bold');
ylabel('D_i (conservation)','FontSize',10,'FontWeight','bold');

% The plots show a histogram of positional conservation values, and a bar
% graph of conservation values for each position.

%% Step 4. SCA calculations

% SCA v5.0 takes a multiple sequence alignment of a protein family as input
% and returns two main outputs: (1) a positional correlation matrix (Cp),
% which quantitatively indicates the correlated evolution of all pairs of
% positions in the alignment, and (2) a sequence correlation matrix (Cs),
% which indicates the pattern of similarity between all pairs of sequences.

[Gsca]=sca5(algn);

%% Step 5. Spectral (or eigenvalue) decomposition

% The function spectral_decomp.m carries out eigenvalue decomposition of
% the Cp matrix and returns a structure with the results for actual and
% randomized alignments. It also makes a plot of the eigenspectra.

[spect]=spectral_decomp(Gsca,100);

% In the MATLAB structure spect, spect.lbdpos and spect.lbdposrnd contain
% the eigenvalues of the actual and randomized alignments, and spect.evpos
% and spect.evposrnd contains the corresponding eigenvectors. For the G
% protein family, it appears that the top ~6 eigenvalues can be considered
% statistically significant given the random matrix expectation (red line).

%%  Step 6. Structure of top eigenmodes

% Here, we examine the structure of the three top eigenmodes here of the Cp
% matrix.

% a 3-D plot of the top three eigenvectors
h_3Dtopmodes=figure; set(h_3Dtopmodes,'Units','normalized','Position',[0 0.7 0.3 0.4],'Name','Top Eigenmodes - 3D'); clf; 
scatter3(spect.evpos(:,1),spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats_1q21);text(spect.evpos(i,1)+.01,spect.evpos(i,2)+.01,spect.evpos(i,3)+.01,ats_1q21(i));end;hold off
az=118;el=12;view(az,el);
xlabel('ev 1','FontSize',12,'FontWeight','b');
ylabel('ev 2','FontSize',12,'FontWeight','b');
zlabel('ev 3','FontSize',12,'FontWeight','b');

% In the case of the G protein family, we find a pattern of positional
% weights on the top modes that shows little evidence of independent
% sectors.  Indeed, the pattern is much like that for the PDZ domain (see
% "Tutoral_PDZ.m"), suggesting a single hierarchically organized sector.
% Accordingly, we will define the sector based on the first eigenmode of
% the Cp matrix.


%% Step 6d. Sector definition

% As in Tutorial_PDZ.m, we fit top eigenvector to the lognormal
% distribution and define the sector by making a cumulative density
% function (cdf) from the fitted distribution and picking a cutoff in the
% tail (at 80%). So, this means we are choosing to represent the residues
% belonging to the top 20% of the cdf.  This cutoff illustrates the basic
% properties of sectors...a few positions in the tail of the distibution of
% eigenvector weights which (as we show further below) comprise distributed
% physically contiguous networks of atoms in the protein structure linking
% in this case the nucleotide binding site to distantly positioned effector
% interaction surfaces.

% figure and probability cutoff
h_secdef=figure;
set(h_secdef,'Units','normalized','Position',[0 1 .25 0.3],'Name','Top Eigenmode'); clf;
p_cutoff=0.8; % the cdf cutoff
secpos = [];

% histogram of the data
binwidth=2*iqr(spect.evpos(:,1))*(numel(spect.evpos(:,1))^-0.33);  % the Freedman-Diaconis rule
nbins=round(range(spect.evpos(:,1))/binwidth);
[yhist,xhist]=hist(spect.evpos(:,1),nbins); bar(xhist,yhist/N_pos,'k');hold on;grid on

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
sprintf('%g+',str2num(char(ats_1q21(sec.def))))
sec.cutoff=cutoff_ev;
figure(h_secdef); line([cutoff_ev cutoff_ev],[0 max(yhist)/N_pos],'LineWidth',1,'LineStyle','--','Color','b');
sec.col=2/3;

%% Step 7a. The mapping between positional and sequence correlations - I

% The simple analysis of unweighted sequence correlations in Step 2 clearly
% suggests the existence of strong subfamily architecture in the G protein
% fmaile.  To see this more clearly, we examine the top eigenmodes of the
% SCA-weighted sequence correaltion matrix Cs. 

[spect.evseq,spect.lbdseq]=eigenvect(Gsca.Cs);

% We plot a histogram of eigenvalues of the Cs matrix:
h_eigenspectrum_seq=figure; set(h_eigenspectrum_seq,'Name','Eigenmodes of SCA matrix'); clf; 
subplot(1,2,1)
[yhist,xhist]=hist(spect.lbdseq,N_seq); bar(xhist,yhist,'k'); axis([min(spect.lbdseq) 1.1*max(spect.lbdseq) 0 1.4*max(yhist)]);hold on
xlabel('eigenvalue of Cs (all)','FontSize',10,'FontWeight','bold');
ylabel('number','FontSize',10,'FontWeight','bold');
grid on
% a zoomed version...which is much clearer.  
subplot(1,2,2)
bar(xhist,yhist,'k'); axis([min(spect.lbdseq) 1.1*(spect.lbdseq(2)) 0 20]);hold on
xlabel('eigenvalue of Cs (2-926)','FontSize',10,'FontWeight','bold');
ylabel('number','FontSize',10,'FontWeight','bold');
grid on

% In the plot, the left graph shows a histogram of all eigenvalues of Cs,
% and the right graph shows a zoomed in version that shows eigenvalues 2
% through N_seq (678).  From this, we see that there is a giant first
% eigenmode (at around 1400) followed by just a few other distinct
% eigenmodes.  We will look at a 3-D plot of the top three eigenvectors:

h_3Dtopmodes_seq=figure; set(h_3Dtopmodes_seq,'Units','normalized','Position',[0 0.7 0.4 0.4],'Name','Top Eigenmodes - 3D'); clf; 
scatter3(spect.evseq(:,1),spect.evseq(:,2),spect.evseq(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');hold on
az=112;el=14;view(az,el);
xlabel('evSeq 1','FontSize',12,'FontWeight','b');
ylabel('evSeq 2','FontSize',12,'FontWeight','b');
zlabel('evSeq 3','FontSize',12,'FontWeight','b');

% The result demonstrates that the G protein alignment clearly contains a
% few distinct subfamilies: the mitochondrial translations initiation
% factors separate along evseq(:,1), the Ras/Rho/Rac family members along
% negative weights in evseq(:,2), the ADP ribosylation factors (ARFs) near
% the origin, the translation elongation factors along positive weights in
% evseq(:,2), and the heterotrimeric G proteins along evseq(:,3).

%% Step 7b.  The mapping between positional and sequence correlations - II

% Is the sector related to the subfamily architecture of the G protein
% family?  To make a mapping between the pattern of positional correlations
% described in the top eigenmodes of the Cp matrix (spect.evpos) and the
% pattern of sequence correlations, we carry out the singular value
% decomposition of the projected sequence alignment (pwX), which is
% returned in the MATLAB structure output by sca5.m (here, Gsca.pwX).  See
% Tutorial_S1A.m for a bit more detail on this step.

[U,sv,V]=svd(Gsca.pwX);

% The matrix Pi=U*V' provides a mathematical mapping between the positional
% and sequence correlation matrices and we will apply the Pi matrix to the
% spect.evpos matrix.  The result is spect.evseq_P, the sequence space
% corresponding to the top eigenmodes of the positional correlation matrix:

N_min=min(N_seq,N_pos);
Pi=U(:,1:N_min)*V(:,1:N_min)';
spect.evseq_P=Pi*spect.evpos;

% We plot the three top eigenvectors in spect.evpos and the corresponding
% patterh of sequence divergence in the spect.evseq_P matrix:

h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Sequence Space by Positional Correlations'); clf; 

h_SectSeq(1)=subplot(1,2,1)
scatter3(spect.evpos(:,1),spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats_1q21);text(spect.evpos(i,1)+.05,spect.evpos(i,2)+.05,spect.evpos(i,3)+.05,ats_1q21(i));end;hold off
az=133;el=44;view(az,el);
xlabel('evpos1','FontSize',12,'FontWeight','b');
ylabel('evpos2','FontSize',12,'FontWeight','b');
zlabel('evpos3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2)
scatter3(spect.evseq_P(:,1),spect.evseq_P(:,2),spect.evseq_P(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
az=133;el=44;view(az,el);
xlabel('evseqP1','FontSize',12,'FontWeight','b');
ylabel('evseqP2','FontSize',12,'FontWeight','b');
zlabel('evseqP3','FontSize',12,'FontWeight','b');

% The result is interesting.  Despite considerable subfamily structure in
% the G protein family, the sequence space corresponding to the top
% eigenmodes of the positional correlation matrix has essentailly no
% structure.  That is, the sector is not correlated with any particular
% mode of sequence divergence.  We note that this is in sharp contrast to
% the Hsp70/110 family (see "Tutorial_Hsp70.m") in which there is a clear
% connection between positional correlations and patterns of sequence
% divergence.  We conclude that the sector in the G protein is a global
% property of the family.

%% The "inverse" mapping: from sequence correlations to positional correlations

% In the previous section, we took the strategy of starting with the
% eigenvectors of the SCA positional correlation matrix (spect.evpos) and
% mapping these to the corresponding modes of sequence divergence
% (spect.evseq_P). However, we can also carry out the inverse mapping -
% starting from patterns in the top eigenmodes of the SCA sequence
% correlation matrix (spect.evseq) and mapping the associated pattern of
% positional correlations (spect.evpos_S). Since the sector fails to
% recognize the distinct subfamilies evident in spect.evpos, we can predict
% that the subfamily architecture revealed by spect.evseq should fail to
% identify the sector.  

% To test this, we start with spect.evseq and apply the Pi matrix to
% produce spect.evpos_S, the corresponding space of positional
% correlations:

spect.evpos_S=Pi'*spect.evseq;

h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Positional Correlations by Sequence Space'); clf; 

h_SectSeq(1)=subplot(1,2,1);
scatter3(spect.evseq(:,1),spect.evseq(:,2),spect.evseq(:,3),'ko','SizeData', 50, 'MarkerFaceColor','k');hold on
az=112;el=18;view(az,el);
xlabel('evseq1','FontSize',12,'FontWeight','b');
ylabel('evseq2','FontSize',12,'FontWeight','b');
zlabel('evseq3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2);
not_sec=find(spect.evpos(:,1)<=cutoff_ev);
scatter3(spect.evpos_S(not_sec,1),spect.evpos_S(not_sec,2),spect.evpos_S(not_sec,3),'ko','SizeData', 50, 'MarkerFaceColor','k');hold on;
scatter3(spect.evpos_S(sec.def,1),spect.evpos_S(sec.def,2),spect.evpos_S(sec.def,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold off;
az=112;el=18;view(az,el);
xlabel('evposS1','FontSize',12,'FontWeight','b');
ylabel('evposS2','FontSize',12,'FontWeight','b');
zlabel('evposS3','FontSize',12,'FontWeight','b');

% The left graph shows the same heterogeneous distribution of sequences in
% the top modes of the SCA sequence correlation matrix Cs.  The right graph
% shows the corresponding modes of positional correlation, with sector
% positions colored in blue.  As predicted, there is no evidence that the
% sector positions correspond to any of the sequence heterogeneity present
% in the alignment.  We note that there are a small number of (nonsector)
% positions that emerge along the first component of evpos_S that do
% correspond to the divergence of the mitochondrial translation initation
% factors (along evseq1); the meaning of this is a matter for further
% study.  

%% Step 10. Structural analysis of sectors

% To examine the organization of the sector, we will just work
% interactively with PyMol (www.pymol.org). Load the pdb file called
% DanK-Sse1Model.pdb and set up your PyMol workspace per your preference.
% For example, we set the background to white, show the structure as a
% cartoon, color the two domains differently (e.g. blue-white for the
% residues 3-383 (the ATPase domain) and light-green for residues 384-607
% (the substrate binding domain)), show small molecules such as ATP in
% stick bonds, and represent the sector as spheres with an overlaid surface
% (surface transparency set to 0.4). It helps to set two-sided lighting to
% be on and the image quality as high as reasonable given your
% computational power.

% To create the sector as an object, we will simply write out sector
% residues in a format suitable for cutting and pasting into PyMol.  

sprintf('%g+',str2num(char(ats_1q21(sec.def))))

% Now use the usual object creation commands in PyMol to create the sector
% as a separate object.  For example,
% create sector, 1Q21 and resi <paste sector here>

%% Conclusion

% Overall, we conclude that SCA yields one sector in the G protein family
% that is a homogeneous property of sequences comprising the alignment.
% This is the case despite considerable heterogeneity of sequence
% divergence amongst members of the G protein family.  This conclusion is
% that quite unlike the case of the Hsp70 sector, which maps to one
% functionally relevance mode of sequence divergence that separates
% allosteric and non-allosteric members of that family.

% Structurally, the G protein sector follows the general empirical rule of
% sectors in all proteins: a physically contiguous network that links in
% this case the nucleotide binding site to distantly positioned surfaces
% mediating downstream signaling.

%%

save Outputs/work_G