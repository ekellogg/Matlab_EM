function [out]=sca5(algn,varargin);
% 
% usages: [SCA_struct]=sca5(algn);
%         [SCA_struct]=sca5(algn, varargin(1), varargin(2),varargin(3));,
%           where varargin represents optional parameters as described
%           below
%
% This function carries out SCA for a protein family given a multiple
% sequence alignment (algn) and a number of optional arguments as decribed
% below.
%
% THE APPROACH  
% ********************
% SCA involves the following procedural steps: 

%   (1) conversion of an alignment in ASCII format (M sequences X L
%       positions) to a M x L x 20 three-dimensional binary tensor
%       (out.x3d). x3d(n,i,a) = 1 if sequence n has amino acid a at
%       position i, and is equal to 0 otherwise.

%   (2) multiplication of the alignment by a weight matrix (out.W) such
%       that each amino acid at each position is weighted by a function of
%       its conservation in the entire alignment (the output is called
%       out.wX, the weighted alignment tensor).  This process is a key step
%       in SCA...it assigns the significance of amino acid in each sequence
%       in the alignment by its evolutionary significance.  The default
%       weighting function is the gradient of relative entropy (the measure
%       of conservation used in SCA).

%   (3) projection of the 3D weighted alignment tensor (WX)  to a 2D
%       weighted alignment matrix (pwX).  Ordinarily, the full SCA for a
%       multiple sequence analysis would result in a 4D-tensor (L x L x 20
%       x 20) containing correlations between every pair of amino acids at
%       every pair of positions.  However, the basic principle of SCA is to
%       infer correlations between positions as judged by the combined
%       information from the ensemble of amino acid pairs at each pair of
%       positions. In other words, we wish to dimension reduce the 4D
%       tensor to a 2D matrix of positional correlations (L x L).  
%       
%       How can we do this dimension reduction? An approach is suggested by
%       analyzing each 20 X 20 amino acid correlation matrix corresponding
%       to each pair of amino acids through a mathematical technique called
%       the singular value decomposition (SVD).  In this decomposition, any
%       p x q matrix X can be written as a product of three matrices:
%       X=u*e*v', where e is a diagonal p x q matrix of so-called "singular
%       values" that indicate the quantity of information in X that is
%       caputured and u and v are p x p and q x q orthonormal matrices,
%       repsectively.  The columns of u and v contain "singular vectors"
%       that indicatee the weights for linear combination of rows and
%       columns, respectively, that contributes to each associated singular
%       value.
%       
%       For SCA, it turns out that the correlations in each 20 x 20 amino
%       acid correlation matrix for each pair of positions i and j are
%       essentially completely described by just their top singular value;
%       thus, a matrix of top singular values for each pair of positions is
%       itself a reasonable scheme for reducing the 4D tensor of amino acid
%       correlations to a 2D matrix of positional correlations.  But there
%       is one additional empirical observation that suggests an even more
%       powerful approach: the singular vector associated with the top
%       singular value for any iven position i is basically invariant with
%       regard to the identity of postion j.  That is, position i makes
%       correlations with all positions j using basically the same weighted
%       combination of amino acids.  This is a non-trivial result that
%       argues that the mean top singular vector for position i taken over
%       all other positions j can be used as a "projection vector" for
%       effectively reducing the alignment itself from the 3D tensor wX (M
%       x L x 20) to a 2D matrix (M x L) that we will call pwX, the
%       projected weighted alignment.  

%   (4) From this 2D version of the alignment, it is possible to compute
%       the SCA positional correlation  matrix (Cp) and the SCA sequence
%       correlation matrix (Cs) directly.
% *************************************************
%
% INPUTS
% ***********************
%
% With the exception that the alignment must be the first argument, the
% order of the subsequent optional arguments does not matter.
%
% - algn    : the multiple sequence alignment
%
% - weights : function handle to a weighting function, which by default is
%             DerivEntropy.m. In general, the weighting function must take
%             two variables as arguments: 'f', the frequency of an amino
%             acid and 'q', the backgroung frequency of that amino acid.
%             The function can optionally also take a single additional
%             parameter called 'para', which is passed to the weighting
%             function.  
% - para    : optional parameter for the weighting function, set by default
%             to para=[];
% - nv_flag : To suppress all screen output during execution, include the
%              input argument 'nv'
%
% ***********************
%
% OUTPUTS 
% *********************** 
%
% The output is a MATLAB structure that contains the following fields: algn
% (the M x L orginal input alignment), X3d (the M x L x 20 binary alignment
% tensor), wX (the M x L x 20 weighted alignment tensor), pwX (the M x L
% projected weighted alignment matrix), pm (the L x 20 projection matrix),
% W (the L x 20 weight matrix), weightfn (the weighting function handle),
% para (the value of the optional parameter, Cp (the SCA positional
% correlation matrix), and Cs (the SCA sequence correlation matrix).
%
% Authors: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%          Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 8/2011
%
% Copyright R.Ranganathan 1999-2011
%**************************************************************************


% Default arguments
weights=@DerivEntropy;
para=[];
nv_flag=0;

% Treatment of optional arguments:
if nargin>1
    for i=2:nargin
        if isa(varargin{i-1},'function_handle'); weights=varargin{i-1};end
        if isa(varargin{i-1},'char'); 
            if strcmp(varargin{i-1},'nv')
                nv_flag=1;
            else
                disp('Option unrecognized...ignoring');
            end
        end
        if isa(varargin{i-1},'numeric'); 
            if strcmp(func2str(weights), 'DerivEntropy')
                disp('No parameters allowed for DerivEntropy weighting function..ignoring para');
            else
                para=varargin{i-1};
            end
        end
    end
end

if nv_flag==0
    fprintf('%s\n',' ');
    disp(['Weights: ' func2str(weights)]);
    if isempty(para)
        disp('Para: None');
    else
        disp(['Para: ' num2str(para)]);
    end
end

% Preliminaries
[N_seq,N_pos]=size(algn);

% ********MAIN PROGRAM *********************
%
% Step 1: Conversion of the MSA into a 3D binary tensor X(a,i,s) with
% x(a,i,s)=1 if and only if sequence s has amino-acid a at position i. The
% function binrep.m does this.

X3d=binrep(algn);

% Step 2: Calculation of the weighted 3d tensor (wX, dimensions N_seqs X
% N_pos X N_aa).  Each amino acid at each position is weighted by a
% function of the conservation of that amino acid in the MSA.  The default
% weighting function in SCA is the gradient of Kullback-Leibler (KL)
% relative entropy (see Note 103), where the KL entropy is the probability
% of observing the frequency of each amino acid at each position given the
% background expectation of that amino acid in protein sequences.  This
% weighting is a core principle of SCA - that the relevance of correlations
% between amino acids scales with the conservation of these amino acids in
% a large and well sampled multiple sequence alignment.  The function
% weight_aln.m does this, and returns the weighted alignment tensor and the
% weight matrix.

[wX,W]=weight_aln(X3d);

% Step 3:  Calculation of the 2D projected alignment matrix (pwX,
% dimensions N_seqs X N_pos) from the 3D weighted alignment tensor (wX).
% See the file header of project_aln.m for details, but in short, the
% projection vector for each position in the MSA is calculated from the
% normalized weighted frequencies of observing each amino acid at each
% position.  Again the weighting function is the gradient of relative
% entropy. The projection matrix (pm, dimensions N_pos X N_aa) contains the
% projection vectors for each position.

[pwX,pm]=project_aln(algn,wX,W);

% Step 4:  Calculation of SCA matrices

Cp=abs((pwX'*pwX)/N_seq-mean(pwX)'*mean(pwX));
Cs=abs((pwX*pwX')/N_pos-mean(pwX')'*mean(pwX'));

% construction of output

out.algn=algn;
out.X3d=X3d;
out.wX=wX;
out.pwX=pwX;
out.pm=pm;
out.W=W;
out.weightfn=weights;
out.para=para;
out.Cp=Cp;
out.Cs=Cs;

end