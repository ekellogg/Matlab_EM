function [w,r]=eigenvect(A,k)
% usage : [w,r]=eigenvect(A)
%         [w,r]=eigenvect(A,k)
%
% Returns all (or the only the first k) principal eigenvectors of A, 
% w(:,1), w(:,2), ..., and their corresponding eigenvalues.
%
% To ensure consistent reproducibility of the calculations, the sign of
% each component, which is otherwise arbitrary, is fixed along the
% direction with most weight.
%
% Author: Olivier Rivoire
% 2/2010
%
%**************************************************************************

if nargin<2, k=size(A,1); end

[v,d]=eig(A);
[l,i]=sort(diag(d),'descend');
w=v(:,i(1:k));
r=l(1:k);

for n=1:k
    w(:,n)=sign(mean(w(:,n)))*w(:,n);
end