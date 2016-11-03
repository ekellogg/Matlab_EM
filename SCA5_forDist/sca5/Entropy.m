function [D, Dmat]=Entropy(f,q,a)
% usages: D=Entropy(f,q)
%         D=Entropy(f,q,a)
%         [D, Dmat]=Entropy(f,q)
%         [D, Dmat]=Entropy(f,q,a)
%
% Gibbs-Shannon Relative entropy D(f|q) if a is not specified or a=1.
% Renyi (relative) entropy with parameter a otherwise.
%
% Inputs:
%
% f can be a scalar 0<f<1 or a vector with f>0 and sum(f)=1.
% q should have same format as for f.
%
% Outputs:
%
% If numel(f)=1, then D is the relative entropy for a particular amino acid
% given its background frequency q and Dmat is undefined.  If numel(f)>1
% (typically, 20, with sum(f)=1), then D is the scalar global relative
% entropy of all amino acids given their respective background frequencies,
% and Dmat is a matrix of individual relative entropies for each amino
% acid.

% Author: Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 2/2010
% revised by R. Ranganathan, 9/2010
%
%**************************************************************************

if nargin<3, a=1; end

if numel(f)==1
    if a==1,
        D=0;
        if f>0, D=D+f*log(f/q); end
        if f<1, D=D+(1-f)*log((1-f)/(1-q)); end
    elseif a==0,
        D=1;
    else
        D=(1/(a-1))*log(f^a*q^(1-a)+(1-f)^a*(1-q).^(1-a));
    end
else % assumes that sum(f)=sum(q)=1 and that f(i)>0 => q(i)>0
    if a==1,
       D=0;
       Dmat=zeros(numel(f),1);
       for i=1:numel(f)
           if f(i)>0 
               D=D+f(i)*log(f(i)/q(i)); 
               if f>0, Dmat(i,1)=Dmat(i,1)+f(i)*log(f(i)/q(i)); end
               if f<1, Dmat(i,1)=Dmat(i,1)+(1-f(i))*log((1-f(i))/(1-q(i))); end
           end
       end
    elseif a==0
        D=1;
    else
        S=0;
        Dmat=zeros(numel(f),1);
        for i=1:numel(f)
            S=S+f(i)^a*q(i)^(1-a);
            Dmat(i,1)=(1/(a-1))*log(f(i)^a*q(i)^(1-a)+(1-f(i))^a*(1-q(i)).^(1-a));
        end
        D=(1/(a-1))*log(S);
    end
end
        