function D=DerivEntropy(f,q)
% usage: D=DerivEntropy(f,q)
%
% Derivative of relative entropy dD(f||q)/df
%
% The default weighting function for SCA is the derivative of relative
% entropy.  This is based on the bootstrap approach for defining the SCA
% positional correlation matrix, as defined in the technical note on SCA
% calculations.  @DerivEntropy is used internally in sca_mat.m.  The
% absolute value is taken to keep  the sign of weighted correaltions
% positive.

% Author: Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 2/2010
% revised 3/2011 by RR
%**************************************************************************

D=0; 
if f>0 & f<1, 
    D=abs(log(f.*(1-q)./(q.*(1-f)))); 
end