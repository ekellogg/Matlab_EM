function [varargout]=cons(algn,varargin);

% usage: [D,Dmat]=cons(algn); 
%
% This function computes the conservation of amino acids at positions in a
% multiple sequence alignment (algn) (see tutorials, or Note 103).  In
% addition, the function outputs a global measure of conservation of all
% amino acids at each position.  In SCA, conservation is computed as the
% Kullback-Leibler (K-L) relative entropy of observing a frequency of each
% amino acid at a position given the backgroud expectation of that amino
% acid in all proteins in the non-redundant database.  The function
% Entropy.m makes this calculation. 
%
% Outputs:  
%
% D :    the positional relative entropy including all amino acids
%        and gaps. Inclusion of gaps is important to make the summed
%        frequencies at each position equal to 1.  The backgroud frequency
%        of gaps is computated from the alignment overall, so it is key to
%        not trivially over-represent gaps.  We use a general rule to
%        impose a maximum gap frequency at positions of 20% (see Note xx).
% Dmat : the matrix of relative entropies of each amino acid at each
%        position.  This is computed without considering gaps.
%
% Inputs: The alignment is ASCII format (algn), and optional arguments in
% varargin.  These could be an alternate function to use for calculation of
% conservation (passed as a function handle), and a single optional
% parameter that if specified will be passed to the conservation function.
% The defaults are K-L relative entropy and no parameter.  Note that
% Entropy.m is equipped to also calculate Renyi entropies if a parameter a
% is specified (see Entropy.m for details).
%
%
% Authors: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%          Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 2/2010, revised 8/2011
% copyright Olivier Rivoire and Rama Ranganathan, 2008-2011.
%**************************************************************************

% convect alignment from text file to numeric
algn=lett2num(algn);

% Default arguments
weights=@Entropy;
para=[];
N_aa=20;

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


fprintf('%s\n',' ');
disp(['Weights: ' func2str(weights)]);
if isempty(para)
    disp('Para: None');
else
    disp(['Para: ' num2str(para)]);
end


% Frequencies and background frequencies, accounting or not for gaps:
prob_ref=[.073 .025 .050 .061 .042 .072 .023 .053 .064 .089 .023 .043 .052 .040 .052 .073 .056 .063 .013 .033];  
prob_gaps=sum(algn(:)==0)/numel(algn(:));
prob_ref_gaps=[(1-prob_gaps)*prob_ref prob_gaps];
[N_seq,N_pos]=size(algn);
for a=1:N_aa, freq(:,a)=sum(algn==a)./size(algn,1); end
freq_wgaps=freq;
freq_wgaps(:,21)=sum(algn==0)./size(algn,1);


% Conservation vectors:
Dglo=zeros(1,N_pos);
Dmat=zeros(20,N_pos);
if isempty(para)
    for p=1:N_pos 
            [Dglo(p)]=feval(weights,freq_wgaps(p,:),prob_ref_gaps);
            for k=1:N_aa;[Dmat(k,p)]=feval(weights,freq(p,k),prob_ref(k));end
            varargout{1}(p)=Dglo(p);
            varargout{2}(:,p)=Dmat(:,p);
    end
else
    for p=1:N_pos,
            [Dglo(p)]=feval(weights,freq_wgaps(p,:),prob_ref_gaps,para);
            for k=1:N_aa;[Dmat(k,p)]=feval(weights,freq(p,k),prob_ref(k),para);end
            varargout{1}(p)=Dglo(p);
            varargout{2}(:,p)=Dmat(:,p);
    end
end

