% SCA 5.0 Tutorial 2- Analysis of an alignment of the S1A family of Serine 
% Proteases.

% Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)
% Olivier Rivoire  (olivier.rivoire@ujf-grenoble.fr)

% - Aug 2011

% This tutorial provides a practical walkthrough of the usage of the SCA
% Toolbox 5.0 using a multiple sequence alignment of the S1A serine
% protease family as an example. The tutorial is in "cell" mode in MATLAB
% in which the commands corresponding to each section (or cell) can be
% executed by selecting the cell and clicking on the MATLAB command for
% "evaluate current cell". This tutorial builds on several concepts
% introduced in 'Tutorial_pdz'. As such, we suggest that the pdz tutorial
% be examined carefully before working through this protein family.

addpath sca5
clear; close all

%% Step 1. Alignment loading 

% The toolbox contains a sample multiple sequence alignment of the S1A
% domain family (comprising 1206 sequences) in .free format, which to be
% consistent with Halabi et al. is already truncated to the structure in
% PDB 3TGI (rat trypsin). Load it using 'get_seqs'. Other file formats
% (e.g. FASTA or MSF) can be loaded using standard functions in the
% Bioinformatics Toolbox.

algn=get_seqs('Inputs/al_S1A_1388.free');
[N_seq,N_pos]=size(algn);

% We make a residue numbering list that maps alignment positions to the
% numbering system of pdb 3TGI. See the "Tutorial_pdz.m" for details on
% this step.

pdb_id='3TGI'; chain='E';
pdb=pdbread(['Inputs/' pdb_id '.pdb']);
[strseqnum,ats,best_align]=MSAsearch(pdb,chain,algn);

% 'strseqnum' identifies the alignment sequence that most closely matches to
% the query molecule (here chainID 'E' of pdb_3TGI), 'ats' represents the
% alignment positions in pdb_3tgi numbering, and 'best_align' shows the
% pairwise alignment between the pdb sequence and the top-hit in the
% alignment. MSAsearch can also be used to truncate the alignment to 
% positions contained in the query molecule (see file header for usage). 
% Here, the alignment is pre-trimmed to the 3tgi structure and we are 
% simply using MSAsearch to make a residue numbering list that corresponds 
% to the alignment positions.

% We load the variables containing the position labels and annotations for
% each sequence in the alignment, for use below.

load Inputs/annot_S1A_1388.mat
load Inputs/labels_sprot_3tgi_223.mat

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

% The histogram shows a reasonably narrow distribution with a mean pairwise
% identity between sequences of about 25% and a range of 10 to 40%,
% suggesting that most sequences are about equally dissimilar from other.
% We can further examine this assertion by direct visualization of the
% sequence similarity (S) matrix:

figure(h_seqsim); 
subplot(1,2,2); imshow(S,[0 1],'InitialMagnification','fit'); colormap(jet); colorbar;
title('SeqID', 'FontSize',12,'FontWeight','bold');

% The similarly matrix basically recapitulates what the histogram tell us;
% there are a few small clades of more related PDZ sequences, but in
% general, the alignment is comprised of a diverse and largely
% homogeneously diverged ensemble of sequences.  

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

[S1Asca]=sca5(algn);
%% Step 5. Spectral (or eigenvalue) decomposition

% The function spectral_decomp.m carries out eigenvalue decomposition of
% the Cp matrix and returns a structure with the results for actual and
% randomized alignments. It also makes a plot of the eigenspectra.

[spect]=spectral_decomp(S1Asca,100);

% In the MATLAB structure spect, spect.lbdpos and spect.lbdposrnd contain the
% eigenvalues of the actual and randomized alignments, and spect.evpos and
% spect.evposrnd contains the corresponding eigenvectors. For the S1A family,
% it appears that just the top 7-10 eigenvalues are statistically
% significant.  These top eigenmodes can be examined for any patterns of
% positional correlations.  

%%  Step 6. Structure of top eigenmodes

% Sectors are empirically defined by examining the pattern of positional
% contributions to the top few eigenvectors. We will just examine the
% structure of the three top eigenmodes here.

% a 3-D plot of the top three eigenvectors
h_3Dtopmodes=figure; set(h_3Dtopmodes,'Units','normalized','Position',[0 0.7 0.3 0.4],'Name','Top Eigenmodes - 3D'); clf; 
scatter3(spect.evpos(:,1),spect.evpos(:,2),spect.evpos(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(spect.evpos(i,1)+.01,spect.evpos(i,2)+.01,spect.evpos(i,3)+.01,ats(i));end;hold off
az=79;el=10;view(az,el);
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

% In the case of the serine proteases, we find a pattern of positional
% weights on the top modes that suggest the existence of multiple at least
% quasi-independent sectors.  For example, a scatter plot of eigenvectors 2
% and 3 shows the emergence of three separate groups of positions along
% different combinations of these two eigenvectors.

% The suggestion of quasi-independent sectors is further explored in the
% next section.  We note that the lower (but still significant) eigenmodes
% contain further more weakly organized groupings of positions.  The
% meaning and interpretation of these weaker clusters is a matter for
% further analysis and is not considered here.

%% Step 7.  Independent Component Analysis (ICA)

% As discussed in Halabi et al.(Cell (2009), 138, 774-786) and in Smock et
% al. (Mol. Sys. Biol. 6, 414), eigenvectors are not expected to directly
% represent statistically independent groups of correlated positions
% (sectors).  In general, if independent sectors exist for a particular
% protein family, they will correspond to groups of positions emerging
% along combinations of eigenvectors (see Step 6 and Halabi et al. for
% example).  The reason is due to the fact that decorrelation of positions
% (by diagonalizing the SCA correlation matrix - the essence of eigenvalue
% decomposition) is a weaker criterion than achieving statistical
% independence (which requires absence of not only pairwise correlations,
% but lack of any higher order statistical couplings).  In other words, if
% the non-independence of a set of variables is not completely captured in
% just their pairwise correlations, then just the linear combination of
% these variables indicated by eigenvectors of the correlation matrix
% cannot be assumed to produce statistically independent transformed
% variables. Many discussions of this problem are available in the
% literature (for example, see the book by Hyvarinen et al. (Independent
% Component Analysis).

% Independent component analysis (ICA) - an extension of spectral
% decomposition - is a heuristic method designed to transform the k top
% eigenmodes of a correlation matrix into k maximally independent
% components through an iterative optimization process.  In principle, this
% process should help to better define independent sectors (if they exist
% in the protein family under study) as groups of positions now projecting
% specifically along the transformed axes - the so-called independent
% components. It is worth noting that ICA assumes the existence of
% quasi-independent components.  In this regard, protein families in which
% multiple independent sectors are not evident (for example, due to the
% case of just one sector) should be associated with failure of ICA to
% clearly separate groups of positions along independent components.

% For the S1A serine proteases, we use ICA to transform the significant top
% eigenmodes of the SCA matrix.  In the SCA toolbox, ICA is computed using
% the simple implementation of the algorithm by Bell and Sejnowski in the
% function basic_ica.m. The variable kmax sets the number of top eigenmodes
% retains for analysis below and basic_ica.m returns the rotation matrix
% (W) used to transform the eigenvectors into independent components (IC).
% We will make the variable containing the independent components "ic_P" to
% indicate that these are independent components of the positional
% correlation matrix Cp.

kmax =8;
learnrate=.0001; iterations=20000;
[W,changes_s]=basic_ica(spect.evpos(:,1:kmax)',learnrate,iterations); 
ic_P=(W*spect.evpos(:,1:kmax)')'; 

% We plot the first three ICs.  However, note that upon ICA transformation,
% the order of IC's does not have obvious meaning with regard to
% significance.  This is in contrast to eigenvectors, which are ordered by
% their corresponding eigenvalue according to quantity of variance
% captured.  Nevertheless, the groupings of amino acid positions in the top
% eigenmodes is found in the top independent components.  For example, the
% three main quasi-independent sectors in the serine protease family are
% found as the top three independent components.  The lower independent
% components may also contain important functional interpretations, but for
% the purposes of this tutorial, we do not consider this here.

h_ICA3D=figure; set(h_ICA3D,'Units','normalized','Position',[0 0.1 0.25 0.4],'Name','Independent Components - 3D'); clf; 
scatter3(ic_P(:,1),ic_P(:,2),ic_P(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(ic_P(i,1)+.05,ic_P(i,2)+.05,ic_P(i,3)+.05,ats(i));end;hold off
az=125;el=42;view(az,el);
xlabel('IC 1','FontSize',12,'FontWeight','b');
ylabel('IC 2','FontSize',12,'FontWeight','b');
zlabel('IC 3','FontSize',12,'FontWeight','b');

% 2-D plots of the top three eigenvectors
h_ICA2D=figure; set(h_ICA2D,'Units','normalized','Position',[0 1 0.75 0.3],'Name','Top Eigenmodes-2D'); clf; 
subplot(1,3,1);
scatter(ic_P(:,1),ic_P(:,2),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(ic_P(i,1)+.01,ic_P(i,2)+.01,ats(i));end;hold off;grid on
xlabel('ev 1','FontSize',12,'FontWeight','b');ylabel('ev 2','FontSize',12,'FontWeight','b');
subplot(1,3,2);
scatter(ic_P(:,1),ic_P(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(ic_P(i,1)+.01,ic_P(i,3)+.01,ats(i));end;hold off;grid on
xlabel('ev 1','FontSize',12,'FontWeight','b');ylabel('ev 3','FontSize',12,'FontWeight','b');
subplot(1,3,3);
scatter(ic_P(:,2),ic_P(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(ic_P(i,2)+.01,ic_P(i,3)+.01,ats(i));end;hold off;grid on
xlabel('ev 2','FontSize',12,'FontWeight','b');ylabel('ev 3','FontSize',12,'FontWeight','b');

% Upon ICA transformation, the data indicate three distinct groups of
% evolutionarily correlated sequence positions (sectors) that emerge along
% the top three independent components. The extent of contribution of each
% position to each sector is given by its value along each independent
% component, and below, we use a quantitative but empirical approach to
% defining the positions that comprise each sector.

% The three sectors are largely independent, but as pointed out in Halabi
% et al., note that some positions jointly contribute.  For example,
% positions 213, 214, and 216 display about equal weights for the sectors
% defined by ICs 1 and 3, and postion 184 seems to display significant
% weight for ICs 1 and 2.

%% Step 8. A mapping between positional and sequence correlations

% As discussed in Tutorial_PDZ, SCA version 5 (sca5) has the ability to map
% between the space of position correlations (whose top eigenmodes define
% the sectors) and the space of sequence correlations (whose top eigenmodes
% tell us about the sub-family architecture, if any) in the alignment.
% Here, we use this mapping to study how the amino acid motifs in the
% serine protease sectors are related to functional or phylogenetic
% differences between sequences.

% To make this mapping between sequence correlations and positional
% correlations, we carry out the singular value decomposition of the
% project sequence alignment (pwX), which is returned in the MATLAB
% structure output by sca5.m (here, S1Asca.pwX).

[U,sv,V]=svd(S1Asca.pwX);

% The important concept is that the matrix Pi=U*V' provides a mathematical
% mapping between the positional and sequence correlation matrices.  We
% used ICA to rotate the top eigenmodes of the Cp matrix to a more
% maximally independent representation, and so we will apply the Pi matrix
% to the IC's of the Cp matrix (contained in ic_P).  The result is U_p, the
% sequence space corresponding to the independent components of the
% positional correlation matrix:

N_min=min(N_seq,N_pos);
Pi=U(:,1:N_min)*V(:,1:N_min)';
U_p=Pi*ic_P;

% Thus, if an eigenmode of the positional correlation matrix (a columns in
% spect.evpos) describes a sector (a group of coevolving amino acid
% positions), then the corresponding column of U_p will contain the
% eigenmode of sequence correlations that will reveal the pattern of
% sequence divergence (if any) that underlies the sector definition.  We
% can make this mapping:

h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Sequence Space by Positional Correlations'); clf; 

h_SectSeq(1)=subplot(1,2,1)
scatter3(ic_P(:,1),ic_P(:,2),ic_P(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(ic_P(i,1)+.01,ic_P(i,2)+.01,ic_P(i,3)+.01,ats(i));end;hold off
az=125;el=42;view(az,el);
xlabel('ev 1','FontSize',12,'FontWeight','b');
ylabel('ev 2','FontSize',12,'FontWeight','b');
zlabel('ev 3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2)
scatter3(U_p(:,1),U_p(:,2),U_p(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
az=125;el=42;view(az,el);
xlabel('Seq 1','FontSize',12,'FontWeight','b');
ylabel('Seq 2','FontSize',12,'FontWeight','b');
zlabel('Seq 3','FontSize',12,'FontWeight','b');

% For the S1A serine proteases, this mapping shows that the sectors do not
% correspond to obvious patterns of inhomogeneities in the sequence space;
% indeed, the top ICA-rotated top eigenmodes of the sequence space are
% largely non-descript. In other words, the sectors do not correspond to
% phylogenetically distinct subfamilies.  We define the sectors below based
% on the three top independent components (ev.ica 1-3) and then we will
% come back and show that there is a meanginful functional partitioning of
% sequences in this family that indicates independent control of different
% aspects of fitness in the serine proteases.


%% Step 9. Sector definition

% For the S1A serine proteases domain, we will define three sectors,
% defined by each of the first three ICs.  The distribution of positional
% weights along independent components comprises a heavy-tailed
% distribution which is well-fit by the t-distribution. Again, as explained
% in Tutorial_PDZ.m, no deep mechanistic basis for this distribution is
% implied here. Instead, we use this empirical fitting simply to provide a
% logical and systematic basis for representing sector positions. 

% We fit the first three ICs to the t location scale distribution and
% define sector positions.  Here, we use the MATLAB function fitdist to
% create so-called probability distribution objects for the fits to each IC
% that provide a convenient way to plot and manipulate the fitted
% distributions.

% we plot the histograms of the ICs along with the fitted distributions. To
% identify sectors, we provide two empirical principles based on study of
% many protein families: (1) sectors are generally defined by roughly
% 10-20% of sequence positions that show physical connectivity in the
% structure and functional importance, (2) sectors should be thought of as
% hierarchically organized and so arguments about the biological relevance
% of sectors should not strongly depend on specific details of statistical
% cutoffs used for their identification.

h_ICAfit=figure;
set(h_ICAfit,'Units','normalized','Position',[0 1 1 0.5],'Name','IC distributions'); clf;
clear sec cutoffs

p_cutoff=0.9;
nfit=3;
cutoffs = zeros(nfit,1);

for i=1:nfit
    pd=fitdist(ic_P(:,i),'tlocationscale');
    subplot(1,nfit,i);
    binwidth=2*iqr(ic_P(:,i))*(numel(ic_P(:,i))^-0.33);  % the Freedman-Diaconis rule
    nbins=round(range(ic_P(:,i))/binwidth);
    % here we plot the histogram of IC weights as probability densities in
    % each bin:
    [yhist,xhist]=hist(ic_P(:,i),nbins); bar(xhist,yhist/N_pos,'k');hold on;grid on
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
    if abs(max(ic_P(:,i)))>abs(min(ic_P(:,i)))
        tail(maxpos:end)=cdf_jnk(maxpos:end);
    else
        cdf_jnk=1-cdf_jnk;
        tail(1:maxpos)=cdf_jnk(1:maxpos);
    end
    [~,x_dist_pos]=min(abs(tail-p_cutoff));
    cutoffs(i) = x_dist(x_dist_pos);

    y_lim=get(gca,'ylim');
    line([cutoffs(i) cutoffs(i)],y_lim,'Color','k','LineWidth',1,'LineStyle','--');
    text(cutoffs(i),1.03*(y_lim(2)-y_lim(1)),num2str(cutoffs(i),'%0.2f'),'FontWeight','bold','FontSize',11);
    xlabel(['IC ' num2str(i)],'FontSize',12,'FontWeight','b');
    ylabel('Prob Density','FontSize',12,'FontWeight','b');
    
    % we obtain the indicies of sector positions
    if abs(max(ic_P(:,i)))>abs(min(ic_P(:,i)))
        sec(i).def = find(ic_P(:,i)>cutoffs(i)); 
    else
        sec(i).def = find(ic_P(:,i)<cutoffs(i)); 
    end
end

% if one wants to gradient color...

% sec(1).vec=ic(sec(1).def,1)-cutoffs(1);
% sec(2).vec=ic(sec(2).def,2)-cutoffs(2);
% sec(3).vec=ic(sec(3).def,3)-cutoffs(3);

% set sector colors

sec(1).col=0;
sec(2).col=2/3;
sec(3).col=1/3;


%% Step 10. Structural analysis of sectors
 
% One of the central findings is that sectors are not obviously related to
% current classification of proteins by primary, secondary, and tertiary
% structure. Indeed, sectors comprise positions that can be widely
% distributed in primary structure, that do not follow known patterns of
% contacts due to secondary structure, and do not follow tertiary
% structural motifs and subdomain architectures. Instead, they appear to be
% new features of three-dimensional structure that follow only the rule
% that they are physically connected in at least one physiological state of
% the protein. To see this in context of the serine protease sectors, we
% map sector positions by primary, secondary, and tertiary structure
% representations.

% Organization by the primary/secondary structure:

h_SecStr=figure;
set(h_SecStr,'Units','normalized','Position',[0 0 0.9 0.35],'Name','Sectors: Primary and Secondary Structure');clf
subplot(2,1,1);
show_vect(D_glo,sec,ats);
SPss.Sheet=[1:13]; SPss.Helix=[1:4];
figure(h_SecStr); hsub=subplot(2,1,2);drawSS(pdb,ats,SPss,sec,['r','b','g']);

% Representation on a 3D structure: a script is generated that displays the
% sectors when executed with the software PyMol (www.pymol.org):

filename='Outputs/sectors_sprot';
write_pmlscript(pdb_id,chain,ats,sec,filename);
fid=fopen([filename '.pml'], 'a');
fprintf(fid,'show stick, (resi 15,16,17) & (chain I)\n');
fprintf(fid,'color yellow, (resi 15,16,17) & (chain I)\n');
fprintf(fid,'set stick_radius, 0.4\n');
fprintf(fid,'show stick, (chain B)\n');
fprintf(fid,'util.cbay main and chain I\n');
fprintf(fid, 'color red, sector_1\n');
fprintf(fid, 'color blue, sector_2\n');
fprintf(fid, 'color green, sector_3\n');
% some other cosmetics
fprintf(fid,'set two_sided_lighting, 1\n');
fprintf(fid,'set sphere_quality, 2\n');
fprintf(fid,'set surface_quality, 2\n');
fprintf(fid,'set stick_quality, 10\n');
fprintf(fid,'set cartoon_oval_quality, 10\n');
fprintf(fid,'set cartoon_loop_quality, 6\n');

fprintf(fid,'set_view (\\\n');
fprintf(fid,'-0.655116498,   -0.695889294,   -0.294215441,\\\n');
fprintf(fid,'0.386854321,    0.025535638,   -0.921788275,\\\n');
fprintf(fid,'0.648974240,   -0.717696548,    0.252478868,\\\n');
fprintf(fid,'0.000000000,    0.000000000, -148.182495117,\\\n');
fprintf(fid,'-12.619968414,  -89.585556030,   -3.959166050,\\\n');
fprintf(fid,'122.692764282,  173.672241211,  -20.000000000 )\n');

fclose(fid);

% The program write_pmlscript.m is not particularly sophisticated or ideal
% in taking advantage of features in PyMol to effectively represent
% sectors.  For example, the few positions here that are shared between
% different sectors are arbitrarily colored by order of sector
% specification in the variable sec. A more accurate but manual way to
% examine sectors is simply to write out sector residues in a format
% suitable for cutting and pasting into PyMol and work interactively
% between MATLAB and PyMol to define the sectors as objects in a PyMol
% session.  For example:

sprintf('%g+',str2num(char(ats(sec(1).def))))
sprintf('%g+',str2num(char(ats(sec(2).def))))
sprintf('%g+',str2num(char(ats(sec(3).def))))

% These commands output residue lists for each sector that can be pasted
% into a selection or object definition in PyMol (Note....sprintf will fail
% as written here if residue labels are not numeric!).  To color sectors by
% strength of contribution, we use the excellent web-downloadable Python
% scripts written by Dr. Robert Campbell of Queen's University, Canada
% (http://pldserver1.biochem.queensu.ca/~rlc/) to (1) write out the weights
% for sector residues into the B-factor column (using data2Bfactor.py) and
% (2) to color by a user-specified gradient (using color_b.py).


%% Step 11. Classification of sequences by sectors

% Halabi et al. showed that the three sectors in the serine proteases
% correspond to distinct patterns of phenotypic divergence (catalytic
% specificity, vertebrate or invertebrate origin, and existence of
% catalytic mechanism) in the protein family.  In that paper, this was done
% by looking at the sequence similarity matrices for all annotated
% sequences in the alignment computed for just the positions contributing
% to each sector.  The top eigenmode of these sector-specific sequence
% similarity matrices was shown to separate sequence by pattern of
% phenotypic divergence.  Here, we show that this same result can be
% recapitulated in a different way....by using our mapping between sectors
% and sequences.  Specificlly, we will examine the modes of the sequence
% correlation matrix that correspond to the three sectors.  We showed above
% that this mapping is globally unremarkable; the sequence divergence in
% U_p associated with the three sectors in ICs 1-3 is not obviously
% separating sequence by distinct subfamilies.  Nevertheless, we show here
% that there are interesting and orthogonal phenotypic divergences that
% correspond to each sector.

% To show this, we first identify the subset of sequences for which
% annotation about primary catalytic specificity is available, or which are
% known to be non-enzymatic (haptoglobins). Vertebrate/invertebrate origin
% is annotated for nearly all sequences.

specificities={'chymotrypsin' 'trypsin' 'tryptase' 'kallikrein' 'granzyme'};
col_rr{1}=['b' 'r' 'y' 'm' 'g'];
annotated=[];
for i=1:numel(specificities)
    annotated=[annotated ; strmatch(specificities(i),textdat(:,7),'exact')];
end
annotated=[annotated ; strmatch('haptoglobin',textdat(:,7),'exact')];

% Truncation of the alignment to the annotated subset:

algn_ann=algn(annotated,:);
textdat_ann=textdat(annotated,:);
N_ann=numel(annotated);

ann{1}=specificities; typ(1)=7; col{1}=[2/3 7/8 1/6 1/12 1/3];
ann{2}={'not vertebrate' 'vertebrate'}; typ(2)=1; col{2}=[4/7 0];
col_rr{2}=['k' 'w'];
ann{3}={'haptoglobin'}; typ(3)=7; col{3}=1/2;
col_rr{3}=['c'];
N_cla=numel(ann);
for c=1:N_cla
    for i=1:numel(ann{c})
        class{c}(i).def=strmatch(ann{c}(i),textdat_ann(:,typ(c)),'exact'); 
        class{c}(i).col=col{c}(i);
    end
end

% Here, we define the four sets of positions -  the three sectors and also
% the full sequence.
N_sec=numel(sec);
U_p_ann=U_p(annotated,:);
h_PCASeq=figure;
set(h_PCASeq,'Units','normalized','Position',[0 0 0.7 1],'Name','Sector-based sequence similarities');clf
for s=1:N_sec
    for c=1:N_cla
        clear yhist
        [yhist_tmp,xhist]=hist(U_p_ann(:,s),25);
        subplot(N_sec,N_cla,(s-1)*N_cla+c);
        bar(xhist,yhist_tmp,'k');grid on;
        for i=1:numel(class{c})
            disp([num2str(s) num2str(c) num2str(i)]);
            [yhist(i,:)]=hist(U_p_ann(class{c}(i).def,s),xhist);
        end
        hold on;hbar=bar(xhist,yhist','stacked');hold off;axis auto
        for a=1:numel(class{c})
            set(hbar(a),'FaceColor',col_rr{c}(a));
        end
        if s==1
            axis([-3 3 0 60]);
        elseif s==2
            axis([-3 3 0 50]);
        else
            axis([-3 3 0 80]);
        end
    end
end

subplot(3,3,1);title('IC1 (red) sector, cat spec');
subplot(3,3,2);title('IC1 (red) sector, origin');
subplot(3,3,3);title('IC1 (red) sector, enz mech');
subplot(3,3,4);title('IC2 (blue) sector, cat spec');
subplot(3,3,5);title('IC2 (blue) sector, origin');
subplot(3,3,6);title('IC2 (blue) sector, enz mech');
subplot(3,3,7);title('IC3 (green) sector, cat spec');
subplot(3,3,8);title('IC3 (green) sector, origin');
subplot(3,3,9);title('IC3 (green) sector, enz mech');

% Using a newer approach, this figure is analagous to Figure 6 of Halabi et
% al., Cell 138, 774-786 (2009) showing the pattern of similarities between
% sequences in the top principal component of the sequence distance matrix.
% The three columns color sequences by their annotations according to
% different criteria (1) column1, according to their primary substrate
% specificity (red-trypsin, magenta-kallikrein, yellow - tryptase,
% blue-chymotrypsin, green-granzyme (there are three different
% specificities within grannzymes)), (2) column2, according to their origin
% (vertebrate-white, invertebrate-black), or (3) column3, according to
% presence of a catalytic mechanism (black-enzymes, cyan-haptoglobins). The
% first row shows sequence distances calculated for the IC 1 (red) sector,
% the second row shows distances calculated for the IC 2 (blue) sector, and
% the third row shows distances calculated for the IC 3 (green) sector. The
% data show that the IC 1 (red) sector separates sequences by catalytic
% specificity, the IC 2 (blue) sector separates sequences by vertebrate or
% invertebrate origin, and the IC 3 (green) sector separates sequences by
% existence of catalytic mechanism.  As shown in Halabi et al., the global
% sequence divergence does not separate sequences well by any of these
% categories.

% As another way to see this mapping, we can color the sequences in top
% three modes of the U_p matrix by (1) the specificity, (2) the
% invertebrate or vertebrate origin, and (3) the presence of enzymatic
% mechanism.

h_ColorUp=figure; set(h_ColorUp,'Units','normalized','Position',[0 0.1 1 .4],'Name','Independent Components - 3D'); clf; 

h_ColorUp(1)=subplot(1,3,1)
hdot=scatter3(U_p_ann(:,1),U_p_ann(:,2),U_p_ann(:,3),'ko','SizeData', 50,'LineWidth',0.1, 'MarkerFaceColor','none');
hold on;
az=116;el=30;view(az,el);
        for a=1:numel(class{1})
            h_tmp=scatter3(U_p_ann(class{1}(a).def,1),U_p_ann(class{1}(a).def,2),U_p_ann(class{1}(a).def,3),'ko','SizeData', 50,'MarkerFaceColor',col_rr{1}(a));
        end
        hold off;grid on
xlabel('Up 1','FontSize',12,'FontWeight','b');
ylabel('Up 2','FontSize',12,'FontWeight','b');
zlabel('Up 3','FontSize',12,'FontWeight','b');

h_ColorUp(2)=subplot(1,3,2)
hdot=scatter3(U_p_ann(:,1),U_p_ann(:,2),U_p_ann(:,3),'ko','SizeData', 50, 'MarkerFaceColor','none');
hold on;
az=116;el=30;view(az,el);
        for a=1:numel(class{2})
            h_tmp=scatter3(U_p_ann(class{2}(a).def,1),U_p_ann(class{2}(a).def,2),U_p_ann(class{2}(a).def,3),'ko','SizeData', 50, 'MarkerFaceColor',col_rr{2}(a));
        end
        hold off;grid on
xlabel('Up 1','FontSize',12,'FontWeight','b');
ylabel('Up 2','FontSize',12,'FontWeight','b');
zlabel('Up 3','FontSize',12,'FontWeight','b');

h_ColorUp(3)=subplot(1,3,3)
hdot=scatter3(U_p_ann(:,1),U_p_ann(:,2),U_p_ann(:,3),'ko','SizeData', 50, 'MarkerFaceColor','none');
hold on;
az=116;el=30;view(az,el);
        for a=1:numel(class{3})
            h_tmp=scatter3(U_p_ann(class{3}(a).def,1),U_p_ann(class{3}(a).def,2),U_p_ann(class{3}(a).def,3),'ko','SizeData', 50, 'MarkerFaceColor',col_rr{3}(a));
        end
        hold off;grid on
xlabel('Up 1','FontSize',12,'FontWeight','b');
ylabel('Up 2','FontSize',12,'FontWeight','b');
zlabel('Up 3','FontSize',12,'FontWeight','b');

%% Conclusion

% SCA for the serine protease family reveals three sectors, each of which
% corresponds to a distinct mode of functional variation in the protein
% family, but which do not correspond obviously to distinct phylogenetic
% subfamily structures in the alignment.  Halabi et al. shows that these
% three sectors experimentally correspond to biochemically distinct
% features of the serine protease.

%%

save Outputs/work_sprot