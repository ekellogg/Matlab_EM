function [w,change]=basic_ica(x,r,Niter,B)
% usages: w=basic_ica(x,r,Niter);
%         w=basic_ica(x,r,Niter,B);
%
% Basic ICA algorithm, based on work by Bell & Sejnowski (infomax).
% See "ftp://ftp.cnl.salk.edu/pub/tony/sep96.public"
% for the original code and more explanations.
%
% Inputs:
% x : L*M input matrix where L = # features and M = # samples;
% r : learning rate / relaxation parameter (e.g. r=.0001);
% Niter : number of iterations (e.g. 1000);
% B (optional): batch block size (e.g. B=30) - default: B=M.
%
% Outputs:
% w : unmixing matrix;
% change : record of incremental changes during the iterations.
%
% r, Niter and B should be adjusted to achieve convergence, which should be
% assessed by visualizing 'change' with 'plot(1:Niter,change)';
%
% The input data should preferentially be sphered, i.e., x'*x=1

% Author: Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
% 2010
% revised 8/2011 by RR
%**************************************************************************

[L,M]=size(x);
w=eye(L); % initialization
change=zeros(Niter,1);
if nargin<4, B=M; end

for sweep=1:Niter % iterations
    w_old=w; t=1;
    for t=t:B:t-1+fix(M/B)*B
        u=w*x(:,t:t+B-1);
        w=w+r*(B*eye(L)+(1-2*(1./(1+exp(-u))))*u')*w;
    end
    delta=reshape(w-w_old,1,numel(w(:)));
    change(sweep)=delta*delta';     
end