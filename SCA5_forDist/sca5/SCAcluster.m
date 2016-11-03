function [p,l,sort_order,sorted,h_mat,h_dend]=SCAcluster(C,pos_labels,max_scale,colormap_in);
% usage: [p,l,sort_order,sorted]=SCAcluster(matrix,pos,1.0,jet,1);
%****************SCAcluster.m*****************
%*********************************************
% Author: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%
%
% Two dimensional hierarchical clustering of SCA correlation matrix using
% city-block distance metric and compete linkage. The function takes in the
% correlation matrix, the position labels, a max_scale for linear mapping
% of the color map to correlation values, the colormap, and a flag
% (raw_or_not) that determines whether an unclustered version of the matrix
% is kicked out as well. Returns the distance Crices for positions (pdist
% output, p), the clusters for positions
%(linkage output, l), the sorted indices
% for positions (sort_order), and figures of the clustered matrix, the
% position dendrogram, and if you choose, the unclustered matrix.
%
% 07/2003 - initial
% 01/2005 - modified for SCA2.0 and R14S1
% 08/2008 - modified for SCA3.0/4.0
%
% Copyright R.Ranganathan 1999-2010
%*********************************************
%*********************************************

[x,y]=size(C);
p=pdist(C,'ci');

if x>2
    l=linkage(p,'co');
    [a,b,sort_order,h_dend]=dend_sca(l,pos_labels,1);
else
    sort_order=[1 2];
end
sorted=C(sort_order,sort_order);

if nargin<4
    map='jet';
else
    map=colormap_in;
end

if nargin<3
    map='jet';
    h_mat=figure;imshow(sorted,[0 max(max(C)')],'InitialMagnification','fit');colormap(map);brighten(0);colorbar
else
    map='jet';
    h_mat=figure;imshow(sorted,[0 max_scale],'InitialMagnification','fit');colormap(map);brighten(0);colorbar
end


function [h,T,sort_order,h_dend]=dend_sca(l,labels,pos_flag);
%usage: [h,T,sort_order]=dend_sca(linkage_output, label_matrix, pos_flag)
%example: [a,b,sort_pert_raw]=dend_sca(l_pert_raw,pert_num,0));
%
%****************dend_sca.m************************
%*********************************************
% Author: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%
%modification of dendrogram function in MATLAB.  Suited for SCA matrix clustering
%as given in Suel and Ranganathan, NSB 10, p.59.
%
%The program ttakes in linkage output, l, and labels, and kicks out a dendrogram
%appropriately labelled.  Set pos_flag to be 1 if you want a left-right 
%dendrogram, or 0 if you ant a top-down dendrogram.
%
% Copyright R.Ranganathan 1999-2008
%*********************************************
%*********************************************

if pos_flag==1
    h_dend=figure;
    [h,T,sort_order]=dend_sca_rel13(l,0,'colorthreshold',0,'orientation','right');axis ij;
    set(gca,'YTickLabel',labels(sort_order));
else
    h_dend=figure;
    [h,T,sort_order]=dend_sca_rel13(l,0,'colorthreshold',0,'orientation','top');
    set(gca,'XTickLabel',labels(sort_order));
    set(get(gca,'XLabel'),'Rotation',90.0)
end

% ************internal function***********************
% ****************************************************

function [h,T,v] = dend_sca_rel13(Z,varargin)
%DENDROGRAM Generate dendrogram plot.
%   DENDROGRAM(Z) generates a dendrogram plot of the hierarchical
%   binary cluster tree Z.  Z is an (M-1)-by-3 matrix, generated
%   by the LINKAGE function, where M is the number of objects in the
%   original dataset.
%
%   *********************************************************
%   Modified from the MATLAB version to custom display axis labels by RR.
%
% Copyright R.Ranganathan 1999-2008

m = size(Z,1)+1;
if nargin < 2
    p = 30;
end

if nargin == 2
    p = varargin{1};
end

orientation = 'd';
horz = false;
color = false;
threshold = 0.7 * max(Z(:,3));

if nargin > 2
    if isnumeric(varargin{1})
        p = varargin{1};
        offset = 1;
    else
        p = 30;
        offset = 0;
    end
    
    if rem(nargin - offset,2)== 0
        error('Incorrect number of arguments to DENDROGRAM.');
    end
    okargs = strvcat('orientation','colorthreshold');
    for j=(1 + offset):2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error(sprintf('Unknown parameter name:  %s.',pname));
        elseif length(k)>1
            error(sprintf('Ambiguous parameter name:  %s.',pname));
        else
            switch(k)
            case 1  % orientation
                if ~isempty(pval) & ischar(pval)
                    orientation = lower(pval(1));
                end
                if ~ismember(orientation,{'t','b','d','r','l'})
                    orientation = 'd';
                    warning('Unknown orientation specified, using ''top''.');
                end
                if ismember(orientation,{'r','l'})
                    horz = true;
                end
            case 2  % colorthreshold
                color = true;
                if ischar(pval) 
                    if ~strmatch(lower(pval),'default')
                        warning('Unknown threshold specified, using default');
                    end
                end
                if isnumeric(pval)
                    threshold = pval;
                end
            end
        end
    end
end  
Z = transz(Z); % convert from m+k indexing to min(i,j) indexing.
T = (1:m)';

% if there are more than p node, dendrogram looks crowded, the following code
% will make the last p link nodes as the leaf node.
if (m > p) & (p ~= 0)
    
    Y = Z((m-p+1):end,:);
    
    R = Y(:,1:2);
    R = unique(R(:));
    Rlp = R(R<=p);
    Rgp = R(R>p);
    W(Rlp) = Rlp;
    W(Rgp) = setdiff(1:p, Rlp);
    W = W';
    T(R) = W(R);
    
    % computer all the leaf that each node (in the last 30 row) has
    for i = 1:p
        c = R(i);
        T = clusternum(Z,T,W(c),c,m-p+1,0); % assign to it's leaves.
    end
    
    
    Y(:,1) = W(Y(:,1));
    Y(:,2) = W(Y(:,2));
    Z = Y;
    
    m = p; % reset the number of node to be 30 (row number = 29).
end

A = zeros(4,m-1);
B = A;
n = m;
X = 1:n;
Y = zeros(n,1);
r = Y;

% arrange Z into W so that there will be no crossing in the dendrogram.
W = zeros(size(Z));
W(1,:) = Z(1,:);

nsw = zeros(n,1); rsw = nsw;
nsw(Z(1,1:2)) = 1; rsw(1) = 1;
k = 2; s = 2;

while (k < n)
    i = s;
    while rsw(i) | ~any(nsw(Z(i,1:2)))
        if rsw(i) & i == s, s = s+1; end
        i = i+1;
    end
    
    W(k,:) = Z(i,:);
    nsw(Z(i,1:2)) = 1;
    rsw(i) = 1;
    if s == i, s = s+1; end
    k = k+1;
end

g = 1;
for k = 1:m-1 % initialize X
    i = W(k,1);
    if ~r(i),
        X(i) = g;
        g = g+1;
        r(i) = 1;
    end
    i = W(k,2);
    if ~r(i),
        X(i) = g;
        g = g+1;
        r(i) = 1;
    end
end
[u,v]=sort(X);
% v is the third output value (PERM)
label = num2str(v');

% set up the color

numColors = 1;theGroups = 1;
groups = 0;
cmap = [0 0 1];

if color
    groups = sum(Z(:,3)< threshold);
    if groups > 1 & groups < (m-1)
        theGroups = zeros(m-1,1);
        numColors = 0;
        for count = groups:-1:1
            if (theGroups(count) == 0)
                P = zeros(m-1,1);
                P(count) = 1;
                P = colorcluster(Z,P,Z(count,1),count);
                P = colorcluster(Z,P,Z(count,2),count);
                numColors = numColors + 1;
                theGroups(logical(P)) = numColors;
            end
        end 
        cmap = hsv(numColors);
        cmap(end+1,:) = [0 0 0]; 
    else
        groups = 1;
    end
    
end  


if  isempty(get(0,'CurrentFigure')) | ishold
    figure;
    set(gcf,'Position', [50, 50, 800, 500]);
else
    newplot;
end

col = zeros(m-1,3);
h = zeros(m-1,1);

for n = 1:(m-1)
    i = Z(n,1); j = Z(n,2); w = Z(n,3);
    A(:,n) = [X(i) X(i) X(j) X(j)]';
    B(:,n) = [Y(i) w w Y(j)]';
    X(i) = (X(i)+X(j))/2; Y(i)  = w;
    if n <= groups
        col(n,:) = cmap(theGroups(n),:);
    else
        col(n,:) = cmap(end,:);
    end
end


ymin = min(Z(:,3));
ymax = max(Z(:,3));
margin = (ymax - ymin) * 0.05;
n = size(label,1);

if(~horz)
    for count = 1:(m-1)
        h(count) = line(A(:,count),B(:,count),'color',col(count,:));
    end
    lims = [0 m+1 max(0,ymin-margin) (ymax+margin)];
    set(gca, 'Xlim', [.5 ,(n +.5)], 'XTick', 1:n, 'XTickLabel', label, ...
        'Box', 'off');
    mask = logical([0 0 1 1]); 
    if strcmp(orientation,'b')
        set(gca,'XAxisLocation','top','Ydir','reverse');
    end 
else
    for count = 1:(m-1)
        h(count) = line(B(:,count),A(:,count),'color',col(count,:));
    end
    lims = [max(0,ymin-margin) (ymax+margin) 0 m+1 ];
    set(gca, 'Ylim', [.5 ,(n +.5)], 'YTick', 1:n, 'YTickLabel', label, ...
        'Box', 'off');
    mask = logical([1 1 0 0]);
    if strcmp(orientation, 'l')
        set(gca,'YAxisLocation','right','Xdir','reverse');
    end
end

if margin==0
    if ymax~=0
        lims(mask) = ymax * [0 1.25];
    else
        lims(mask) = [0 1];
    end
end
axis(lims);

% ---------------------------------------
function T = clusternum(X, T, c, k, m, d)
% assign leaves under cluster c to c.

d = d+1;
n = m; flag = 0;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = clusternum(X, T, c, k, n,d); % trace back left subtree
        T = clusternum(X, T, c, X(n,2), n,d);
        flag = 1; break;
    end
end

n = size(X,1);
if flag == 0 & d ~= 1 % row m is leaf node.
    T(X(m,1)) = c;
    T(X(m,2)) = c;
end
% ---------------------------------------
function T = colorcluster(X, T, k, m)
% find local clustering

n = m; 
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = colorcluster(X, T, k, n); % trace back left subtree
        T = colorcluster(X, T, X(n,2), n);
        break;
    end
end
T(m) = 1;
% ---------------------------------------
function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.  

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formmed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

m = size(Z,1)+1;

for i = 1:(m-1)
    if Z(i,1) > m
        Z(i,1) = traceback(Z,Z(i,1));
    end
    if Z(i,2) > m
        Z(i,2) = traceback(Z,Z(i,2));
    end
    if Z(i,1) > Z(i,2)
        Z(i,1:2) = Z(i,[2 1]);
    end
end


function a = traceback(Z,b)

m = size(Z,1)+1;

if Z(b-m,1) > m
    a = traceback(Z,Z(b-m,1));
else
    a = Z(b-m,1);
end
if Z(b-m,2) > m
    c = traceback(Z,Z(b-m,2));
else
    c = Z(b-m,2);
end

a = min(a,c);

