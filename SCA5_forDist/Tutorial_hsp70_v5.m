% SCA 5.0 Tutorial 3 - Analysis of an alignment of the Hsp70/110 family of
% chaperones.

% Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)
% Olivier Rivoire  (olivier.rivoire@ujf-grenoble.fr)

% - Aug 2011

% This tutorial presents the analysis for an alignment - the Hsp70/110
% family - where the sequences are distributed into distinct subfamilies.
% It builds on the concepts introduced in 'Tutorial_pdz' and
% 'Tutorial_sprot'.  As such, we suggest that these tutorials be examined
% carefully before working through this protein family.

% This tutorial is based on work on the Hsp70/110 chaperones by R. Smock,
% O. Rivoire, W. Russ, J. Swain, S. Leibler, R. Ranganathan & L. Gierasch,
% published in Smock et al.(Mol. Syst. Biol. 6, 414. We refer to this
% article for interpreting the results of the analysis.

addpath sca5
clear; close all

%% Step 1.  Alignment and other variable loading

% The toolbox contains a multiple sequence alignment of the Hsp70/110
% domain family (comprising 926 sequences) in .fasta format. Load it using
% 'fastaread.m'. Other file formats (e.g. MSF or .free) can be loaded using
% standard functions in the Bioinformatics Toolbox or by using
% "get_seqs.m", respectively.

algn_tmp=fastaread('Inputs/al_hsp70.fasta');
algn=char({algn_tmp.Sequence});
[N_seq,N_pos]=size(algn);

% We load a structural model of the ATP-bound configuration of Hsp70 and
% load a residue numbering list (ats) that maps alignment positions to the
% numbering system of the structural model. See Smock et al. for details
% about the structural model.

pdb_id='DnaK-Sse1Model'; chain=[];
pdb=pdbread(['Inputs/' pdb_id '.pdb']);
load Inputs/ats_hsp70 ats
for i=1:numel(ats);ats_tmp{i}=num2str(ats(i));end;ats=ats_tmp;clear ats_tmp

% The Hsp70/110s are proteins containing two domains....an N-terminal
% ATPase domain and a C-terminal substrate binding domain (SBD).  We define
% the two domains here.

dom(1).def=[strmatch('3',ats,'exact'):strmatch('383',ats,'exact')]'; dom(1).col=2/3;
dom(2).def=[strmatch('384',ats,'exact'):strmatch('607',ats,'exact')]'; dom(2).col=1/6;
for i=1:2, dom(i).vec=ones(N_pos,1); end

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

% The histogram shows a rather broad and multi-peaked distribution with
% groups of sequence showing mean pairwise identities between sequences
% ranging from 20% to greater than 50%.  This pattern already suggests the
% existence of multiple subfamilies. We can further examine this assertion
% by direct visualization of the sequence similarity (S) matrix:

figure(h_seqsim); 
subplot(1,2,2); imshow(S,[0 1],'InitialMagnification','fit'); colormap(jet); colorbar;
title('SeqID', 'FontSize',12,'FontWeight','bold');

% The similarity matrix again indicates the presence of several subfamilies
% of sequences.  Here, we make use of a clear functional distinction
% between members of the Hsp70 and Hsp110 subfamilies to identify a sector
% connected with this function.  In general, the rich subfamily structure
% in this alignment might tell us more about specific sectors related to
% various aspects of functional divergence between subfamilies, but for the
% purpose of this tutorial we focus only on the main aspect of sequence
% divergence.


%% Step 3. Positional Conservation

% We compute the degree of conservation of
% each position by the Kullback-Leibler relative entropy - D(f(a,i)||q(a))
% See Note 109 and Tutorial_pdz.m for further details. 

[D_glo]=cons(algn);

h_D=figure; set(h_D,'Units','normalized','Position',[0 0.6 0.5 0.4],'Name','Positional Conservation');clf
subplot(2,1,1);hist(D_glo,25); grid on;
xlabel('D (conservation)','FontSize',10,'FontWeight','bold'); 
ylabel('number','FontSize',10,'FontWeight','bold');
subplot(2,1,2);bar([1:numel(ats)],D_glo,'k'); grid on;
axis([0 numel(ats)+1 0 4]);
set(gca,'XTick',[1:10:numel(ats)]);
set(gca,'XTickLabel',ats([1:10:numel(ats)]));
xlabel('position (DnaK numbering)','FontSize',10,'FontWeight','bold');
ylabel('D_i (conservation)','FontSize',10,'FontWeight','bold');

% The plots show a histogram of positional conservation values, and a bar
% graph of conservation values for each position.

%% Step 4. SCA calculations

% SCA v5.0 takes a multiple sequence alignment of a protein family as input
% and returns two main outputs: (1) a positional correlation matrix (Cp),
% which quantitatively indicates the correlated evolution of all pairs of
% positions in the alignment, and (2) a sequence correlation matrix (Cs),
% which indicates the pattern of similarity between all pairs of sequences.

[Hsp70sca]=sca5(algn);

%% Step 5. Spectral (or eigenvalue) decomposition

% The function spectral_decomp.m carries out eigenvalue decomposition of
% the Cp matrix and returns a structure with the results for actual and
% randomized alignments. It also makes a plot of the eigenspectra.

[spect]=spectral_decomp(Hsp70sca,100);

% In the MATLAB structure spect, spect.lbdpos and spect.lbdposrnd contain the
% eigenvalues of the actual and randomized alignments, and spect.evpos and
% spect.evposrnd contains the corresponding eigenvectors. For the Hsp70/110
% family, it appears that the top ~20 eigenvalues can be considered
% statistically significant give the random matrix expectation (red line).
% However, it is clear that the top four eigenvalues stand out
% significantly from the rest, and for the purposes of illustrating sector
% identificaton through functional inhomogeneity in the sequence space (the
% main goal of this tutorial), we will just consider the pattern of
% correlations in the top four modes.  

%%  Step 6. Structure of top eigenmodes

% Here, we examine the structure of the three top eigenmodes here of the Cp
% matrix.

% a 3-D plot of the top three eigenvectors
h_3Dtopmodes=figure; set(h_3Dtopmodes,'Units','normalized','Position',[0 0.7 0.3 0.4],'Name','Top Eigenmodes - 3D'); clf; 
scatter3(spect.evpos(:,1),spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,1)+.01,spect.evpos(i,2)+.01,spect.evpos(i,3)+.01,ats(i));end;hold off
az=118;el=12;view(az,el);
xlabel('ev 1','FontSize',12,'FontWeight','b');
ylabel('ev 2','FontSize',12,'FontWeight','b');
zlabel('ev 3','FontSize',12,'FontWeight','b');

% In the case of the Hsp70/110 family, we find a pattern of positional
% weights on the top modes that suggest the existence of multiple
% quasi-independent sectors.  For example, the scatter plot of eigenvectors
% 1, 2, and 3 shows the emergence of separate groups of positions along
% different combinations of these eigenvectors. In the the Tutorial_S1A.m,
% we showed a similar result with regard to positional correlations, but in
% this case, the corresponding pattern of sequence correlations showed
% little evidence of subfamily architecture.  Thus, we proceeded with
% sector defintion through ICA-based rotation of the top eigenmodes of the
% Cp matrix, and subsequent examination of the functional interpretation of
% the sectors.

% However, the inhomogeneity of sequence relationships in the Hsp70/110
% alignment and the fact that the alignment is known to contain at least
% two functionally distinct sub-families - the allosteric Hsp70 proteins,
% and the non-allosteric Hsp110s - suggests a targeted approach to sector
% identification that is guided by the sub-family architecture. Below, we
% will examine the eigenvalue and ICA decomposition of the Cs matrix
% (instead of the Cp matrix), and will identify a sector corresponding to
% the mode of sequence divergence that separates the Hsp70 and Hsp110
% families.  As we show in Smock et al, MSB 6,414, this sector corresponds to the
% allosteric mechanism in Hsp70.

%% Step 7.  Independent Component Analysis (ICA)

% As discussed in the Tutorial_S1A (and in Halabi et al.(Cell (2009), 138,
% 774-786) and Smock et al. (Mol. Sys. Biol. 6, 414)), eigenvectors are not
% expected to directly represent statistically independent groups of
% correlated positions (sectors).  In the general case, if independent
% sectors exist, they will correspond to groups of positions emerging along
% combinations of eigenvectors.  The reason is due to the fact that just
% decorrelation of positions by eigenvalue decomposition is not guaranteed
% to achieve statistical independence (which requires absence of not only
% pairwise correlations, but lack of any higher order statistical
% couplings).

% As in Tutorial_S1A, we will use independent component analysis (ICA) - an
% extension of spectral decomposition - to transform the k chosen top
% eigenmodes of a correlation matrix into k maximally independent
% components. ICA is computed using the function basic_ica.m. The variable
% kmax sets the number of top eigenmodes retains for analysis below and
% basic_ica.m returns the rotation matrix (W) used to transform the
% eigenvectors into independent components (IC). We will make the variable
% containing the independent components "ica.icpos" to indicate that these are
% independent components of the positional correlation matrix Cp.  We
% choose to apply ICA to the top four eigenmodes:

kmax=4;
learnrate=.0001; iterations=20000;
[W,changes_s]=basic_ica(spect.evpos(:,1:kmax)',learnrate,iterations); 
ica.icpos=(W*spect.evpos(:,1:kmax)')'; 

% Also, to make a mapping between the pattern of positional correlations
% described in ica.icpos (which defines sectors) and sequence correlations, we
% carry out the singular value decomposition of the projected sequence
% alignment (pwX), which is returned in the MATLAB structure output by
% sca5.m (here, Hsp70sca.pwX).  See Tutorial_S1A.m for a bit more detail on
% this step.

[U,sv,V]=svd(Hsp70sca.pwX);

% The matrix Pi=U*V' provides a mathematical mapping between the positional
% and sequence correlation matrices and we will apply the Pi matrix to the
% ica.icpos matrix.  The result is ica.icseq_P, the sequence space
% corresponding to the independent components of the positional correlation
% matrix:

N_min=min(N_seq,N_pos);
Pi=U(:,1:N_min)*V(:,1:N_min)';
ica.icseq_P=Pi*ica.icpos;

% We plot the three top ICs in ica.icpos and the corresponding patterh of
% sequence divergence in the ica.icseq_P matrix:

h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Sequence Space by Positional Correlations'); clf; 

h_SectSeq(1)=subplot(1,2,1);
scatter3(ica.icpos(ica.icpos(:,2)>=-2,1),ica.icpos(ica.icpos(:,2)>=-2,2),ica.icpos(ica.icpos(:,2)>=-2,3),'ko','SizeData', 50, 'MarkerFaceColor','k');hold on
scatter3(ica.icpos(ica.icpos(:,2)<-2,1),ica.icpos(ica.icpos(:,2)<-2,2),ica.icpos(ica.icpos(:,2)<-2,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
for i=1:numel(ats);text(ica.icpos(i,1)+.01,ica.icpos(i,2)+.01,ica.icpos(i,3)+.01,ats(i));end;hold off
az=106;el=50;view(az,el);axis tight
xlabel('icpos1','FontSize',12,'FontWeight','b');
ylabel('icpos2','FontSize',12,'FontWeight','b');
zlabel('icPpos3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2);
scatter3(ica.icseq_P(ica.icseq_P(:,2)<=1,1),ica.icseq_P(ica.icseq_P(:,2)<=1,2),ica.icseq_P(ica.icseq_P(:,2)<=1,3),'ko','SizeData', 50, 'MarkerFaceColor','k');hold on
scatter3(ica.icseq_P(ica.icseq_P(:,2)>1,1),ica.icseq_P(ica.icseq_P(:,2)>1,2),ica.icseq_P(ica.icseq_P(:,2)>1,3),'ko','SizeData', 50, 'MarkerFaceColor','b');hold off
az=106;el=50;view(az,el);axis tight
xlabel('icseqP1','FontSize',12,'FontWeight','b');
ylabel('icseqP2','FontSize',12,'FontWeight','b');
zlabel('icseqP3','FontSize',12,'FontWeight','b');

% Note that colors used here are illustrative to highlight the positional
% correlations and sequence divergences.  Below, we will more
% systematically define the sector by cutoff from a fitted distribution.

% The result shows the main conclusion of sector mapping the Hsp70/110
% family reported in Smock et al. MSB 6, 414. The second mode of ica.icseq_P (see
% the graph at right) shows a clear separation of two groups of sequences,
% a major class near the origin which corresponds to the allosteric
% Hsp70-like sequences, and a small subgroup with positive weights that
% corresponds to the non-allosteric Hsp110-like sequences.  Interestingly,
% the corresponding mode in the ica.icpos matrix (left graph) shows a clear
% separation of a subset of positions with negative weights (the sign if
% independent components is arbitrary) that represent a sector that should
% underlie the Hsp70-110 divergence.  That is, this sector is a hypothesis
% for the positions underlying allostery in the Hsp70 subfamily of
% sequences.

%% Step 8. Sector definition

% We will define the allosteric sector by the distribution along component
% 2 of ica.icpos. To do this, we will fit this IC the t location scale
% distribution and define sector positions by a chosen cutoff, here the top
% 10% of the fitted cumulative density function along the second mode of
% ica.icpos.

h_ICAfit=figure;
set(h_ICAfit,'Units','normalized','Position',[0 1 .4 0.3],'Name','IC distributions'); clf;
clear sec cutoffs

p_cutoff=0.1;
nfit=1;
cutoffs = zeros(nfit,1);

pd=fitdist(ica.icpos(:,2),'tlocationscale');
binwidth=2*iqr(ica.icpos(:,2))*(numel(ica.icpos(:,2))^-0.33);  % the Freedman-Diaconis rule
nbins=round(range(ica.icpos(:,2))/binwidth);
% here we plot the histogram of IC weights as probability densities in
% each bin:
[yhist,xhist]=hist(ica.icpos(:,2),nbins); bar(xhist,yhist/N_pos,'k');hold on;grid on
% we plot the fitted distribution:
x_dist=[min(xhist):(max(xhist)-min(xhist))/100:max(xhist)];
area_hist=N_pos*(xhist(2)-xhist(1)); % for proper scaling of the pdf
pdf_jnk=pdf(pd,x_dist);
scaled_pdf=area_hist.*pdf_jnk;
plot(x_dist,scaled_pdf./N_pos,'r-','LineWidth',1.5);

cdf_jnk=cdf(pd,x_dist);
% here, we identify the direction of the tail (the sign of independent
% components is arbitrary), the cutoff, and the sector positions based
% on the fitted cdf:
[~,maxpos]=max(pdf_jnk);
tail=zeros(1,numel(pdf_jnk));
if abs(max(ica.icpos(:,2)))>abs(min(ica.icpos(:,2)))
    tail(maxpos:end)=cdf_jnk(maxpos:end);
else
    cdf_jnk=1-cdf_jnk;
    tail(1:maxpos)=cdf_jnk(1:maxpos);
end
[~,x_dist_pos]=min(abs(tail-p_cutoff));
cutoffs = x_dist(x_dist_pos);

y_lim=get(gca,'ylim');
line([cutoffs cutoffs],y_lim,'Color','k','LineWidth',1,'LineStyle','--');
text(cutoffs,1.03*(y_lim(2)-y_lim(1)),num2str(cutoffs,'%0.2f'),'FontWeight','bold','FontSize',11);
xlabel(['IC2'],'FontSize',12,'FontWeight','b');
ylabel('Prob Density','FontSize',12,'FontWeight','b');

% we obtain the indicies of sector positions
if abs(max(ica.icpos(:,2)))>abs(min(ica.icpos(:,2)))
    sec.def = find(ica.icpos(:,2)>cutoffs);
else
    sec.def = find(ica.icpos(:,2)<cutoffs);
end

sec.col=2/3;

%% Step 9. Structural analysis of sectors
 
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

sprintf('%g+',str2num(char(ats(sec.def))))

% Now use the usual object creation commands in PyMol to create the sector
% as a separate object.  For example,
% create sector, DnaK-SseModel and resi <paste sector here>

% To color sectors by strength of contribution, we use the excellent
% web-downloadable Python scripts written by Dr. Robert Campbell of Queen's
% University, Canada (http://pldserver1.biochem.queensu.ca/~rlc/) to (1)
% write out the weights for sector residues into the B-factor column (using
% data2Bfactor.py) and (2) to color by a user-specified gradient (using
% color_b.py).

%% An alternative, basically equivalent approach.

% SCA5 returns not only the SCA positional correlation matrix Cp, but also
% the SCA sequence correlation matrix Cs.  In principle, given the results
% above, we should be able to identify the same sector through an "inverse"
% approach from that presented above in which we use spectral decomposition
% and ICA on the Cs matrix to (1) identify the independent component of
% sequence divergence that represents the distinction of Hsp70-like
% sequences and Hsp110-like sequences and then (2) use the Pi projection
% matrix to identify the corresponding component of positional correlation
% that underlies this mode of sequence divergence.  

% To show this, we first carry out an eigenvalue and ICA decomposition of
% Cs matrix:

% eigenvalues and eigenvectors of the Cs matrix, and a plot of the
% eigenspectrum.
[spect.evseq,spect.lbdseq]=eigenvect(Hsp70sca.Cs);
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
% through N_seq (926).  From this, we see that there is a giant first
% eigenmode (at around 3000) followed by just a few other distinct
% eigenmodes.  We will look at a 3-D plot of the top three eigenvectors:

% empirical definition of the different subfamilies evident in the sequence
% space
fam1=find(spect.evseq(:,1)<.025);
fam2=find(spect.evseq(:,2)<-.05);
fam3=find(spect.evseq(:,3)<-0&spect.evseq(:,2)>-.02);
fam4=find(spect.evseq(:,1)>=.025&spect.evseq(:,2)>=-.05&(spect.evseq(:,3)>=0|spect.evseq(:,2)<=-.02));

h_3Dtopmodes_seq=figure; set(h_3Dtopmodes_seq,'Units','normalized','Position',[0 0.7 0.8 0.4],'Name','Top Eigenmodes - 3D'); clf; 
subplot(1,2,1);
scatter3(spect.evseq(fam4,1),spect.evseq(fam4,2),spect.evseq(fam4,3),'ko','SizeData', 50, 'MarkerFaceColor','g');hold on
scatter3(spect.evseq(fam1,1),spect.evseq(fam1,2),spect.evseq(fam1,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
scatter3(spect.evseq(fam2,1),spect.evseq(fam2,2),spect.evseq(fam2,3),'ko','SizeData', 50, 'MarkerFaceColor','r');
scatter3(spect.evseq(fam3,1),spect.evseq(fam3,2),spect.evseq(fam3,3),'ko','SizeData', 50, 'MarkerFaceColor','y');
az=112;el=14;view(az,el);
xlabel('evSeq 1','FontSize',12,'FontWeight','b');
ylabel('evSeq 2','FontSize',12,'FontWeight','b');
zlabel('evSeq 3','FontSize',12,'FontWeight','b');

% These results show that there are four reasonably obvious subfamilies of
% seuqences in the Hsp70/110 alignment.  The red, green, and yellow all
% correspond to Hsp70-like sequences, and the blue corresponds to
% Hsp110-like sequences.  Though clearly separated, the top eigenmodes do
% not provide an optimal representation of these sequence subfamilies; for
% example, some groups diverge from each other on combinations of the top
% eigenmodes.  This is for the same reasons as explained previously for the
% case as for eigenmodes of the positional correlation matrix: in general,
% decorrelation (the essence of eigenvalue decomposotion) does not
% guarantee statistical independence of the eigenmodes modes.  To achieve
% this, more advanced methods such as independent component analysis (ICA)
% can be used, which linearly rotates the k top eigenmodes of a correlation
% matrix to define k maximally independent components.

% Here, we use ICA on the top four eigenmodes of the Cs matrix, and plot.
% We will call the resulting independent components, ica.icseq to denote that
% these are independent components of the SCA-weighted sequence correlation
% matrix Cs.

kmax=4;
learnrate=.0001; iterations=20000;
[W,changes_s]=basic_ica(spect.evseq(:,1:kmax)',learnrate,iterations); 
ica.icseq=(W*spect.evseq(:,1:kmax)')'; 

subplot(1,2,2)
scatter3(ica.icseq(fam4,1),ica.icseq(fam4,2),ica.icseq(fam4,3),'ko','SizeData', 50, 'MarkerFaceColor','g');hold on
scatter3(ica.icseq(fam1,1),ica.icseq(fam1,2),ica.icseq(fam1,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
scatter3(ica.icseq(fam2,1),ica.icseq(fam2,2),ica.icseq(fam2,3),'ko','SizeData', 50, 'MarkerFaceColor','r');
scatter3(ica.icseq(fam3,1),ica.icseq(fam3,2),ica.icseq(fam3,3),'ko','SizeData', 50, 'MarkerFaceColor','y');
az=112;el=14;view(az,el);
xlabel('icseq 1','FontSize',12,'FontWeight','b');
ylabel('icseq 2','FontSize',12,'FontWeight','b');
zlabel('icseq 3','FontSize',12,'FontWeight','b');

% As expected, the ICA algorithm now provides a "better" representation of
% the top eigenmodes of the sequence correlation matrix; the four major
% subfamilies are now represented such that the green group is located at
% the origin and each other sub-family (red, blue, and yellow) emerge
% along orthogonal independent components (ICs).  


%% Step 8. A mapping between positional and sequence correlations

% We now map the positional correlations that underlie these sequence
% divergences by using the Pi matrix, which as explained above provides a
% mapping between positional and sequence correlations.  The result of
% applying the Pi matrix to the independent components of the sequence
% correlations (ica.icseq) is ica.icpos_S, the corresponding space of
% positional correlations:

ica.icpos_S=Pi'*ica.icseq;

% Thus, if an independent component of the sequence correlation matrix (a
% column in ica.icseq) describes a divergence of two subfamilies, then the
% corresponding column of ica.icpos_S will contain the corresponding
% component of positional correlations that underlies this sequence
% divergence.  We can make this mapping:

h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Positional Correlations by Sequence Space'); clf; 

h_SectSeq(1)=subplot(1,2,1);
scatter3(ica.icseq(fam4,1),ica.icseq(fam4,2),ica.icseq(fam4,3),'ko','SizeData', 50, 'MarkerFaceColor','k');hold on
scatter3(ica.icseq(fam1,1),ica.icseq(fam1,2),ica.icseq(fam1,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
scatter3(ica.icseq(fam2,1),ica.icseq(fam2,2),ica.icseq(fam2,3),'ko','SizeData', 50, 'MarkerFaceColor','k');
scatter3(ica.icseq(fam3,1),ica.icseq(fam3,2),ica.icseq(fam3,3),'ko','SizeData', 50, 'MarkerFaceColor','k');
az=112;el=18;view(az,el);
xlabel('icseq1','FontSize',12,'FontWeight','b');
ylabel('icseq2','FontSize',12,'FontWeight','b');
zlabel('icseq3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2);
not_sec=find(ica.icpos(:,2)>=cutoffs);
scatter3(ica.icpos_S(not_sec,1),ica.icpos_S(not_sec,2),ica.icpos_S(not_sec,3),'ko','SizeData', 50, 'MarkerFaceColor','k');hold on;
scatter3(ica.icpos_S(sec.def,1),ica.icpos_S(sec.def,2),ica.icpos_S(sec.def,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold off;
az=112;el=18;view(az,el);
xlabel('icposS1','FontSize',12,'FontWeight','b');
ylabel('icposS2','FontSize',12,'FontWeight','b');
zlabel('icposS3','FontSize',12,'FontWeight','b');

% This result reveals a groups of positional correlations that maps in the
% same direction as the sequence divergence that separates the allosteric
% from non-allosteric members of Hsp70/110 family.  Not surprisingly, these
% positions comprise the same sector as identified by analyzing the SCA
% positional correlation matrix (Cp).

%% Step 10. Conclusion

% For the Hsp70/110 family, we find evidence for multiple sectors in the
% context of an inhomogeneous sequence alignment.  Interestingly, one of
% the top modes of sequence variation in this family cleanly separates the
% allosteric Hsp70-like sequneces from the non-allosteric Hsp110-like
% sequences.  This mode corresponds in positional correlation space to a
% sector comprising roughly 20% of total positions.  Structural mapping
% onto a model for the ATP-bound configuration of the E.coli Hsp70 (DnaK)
% shows that this sector comprises a physically contiguous network of amino
% acids linking the ATP-binding site in the nucleotide binding (or ATPase)
% domain to the ligand bidning site in the substrate binding domain through
% the interdomain interface.  As described in Smock et al. MSB 6, 414, this
% sector is consistent with experimental evidence for the allosteric
% mechanism in Hsp70.

%%
save Outputs/work_hsp70