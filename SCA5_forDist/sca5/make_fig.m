function make_fig(v1,v2,sec)
% usage: []=make_fig_col(v1,v2,sec);
%
% 2d scatter plot with optional gradients of color and size for the points.
% 
% v1 and v2 should be 2 vectors of same size: they specify the coordinates 
% of the points.
%
% 'sec' should be a structure with fields:
% - def
% - vec
% - col
% - shp
% - siz
%
% Only the field 'def' is required and any other field may be omitted, 
% in which case it will be assigned a default value. Omissions should
% however be consistent accross all the entries of 'sec'.
%
% The color code is visualized with 'make_fig(0)' : the color at angle
% alpha is coded by alpha/2*PI (alpha in radian).
% For instance: 0->red, 1/3->green, 2/3-> blue.
% Admissible color codes are in [0 1] with 1 coding for white (gray can be 
% obtained by using -1).
% 
% Example:
%{
    v1=rand(100,1)-.5; v2=rand(100,1)-.5;
    sec(1).def=find(v1>0); sec(1).vec=v1; sec(1).col=0; 
    sec(1).shp='s'; sec(1).siz=v2; 
    sec(2).def=find(v1<0); sec(2).vec=-v1; sec(2).col=2/3;
    sec(2).shp='d';  sec(2).siz=v2;
    figure(1); clf; make_fig(v1,v2,sec);
%}
%
% Note: the order in 'sec' matters: each group of points is plotted
% successively and some points may be masked by others.

% Author: Olivier Rivoire (orivoire@rockefeller.edu)
% 2/2010
%
%**************************************************************************

%% Colormap

% If the only argument is '0', returns the color map for reference
if nargin==1 & v1==0, 
    pi=3.14159;
    hold on;
    for s=0:.05:1
        for a=0:.01:1
            [b,g,r]=hsv2rgb([a,s,1]);
            plot(s*cos(2*pi*a),s*sin(2*pi*a),'o','MarkerSize',8,...
                'MarkerFaceColor',[b g r],'MarkerEdgeColor',[b g r]);
        end
    end
    axis([-1.1 1.1 -1.1 1.1]);
    return;
end

%% Completion of optional fields

% Number of points and number of sectors
N_pts=numel(v1); N_sec=numel(sec);

% Default values for non-specified fields
if ~isfield(sec,'vec')
    for i=1:N_sec, sec(i).vec=ones(N_pts,1); end
end
if ~isfield(sec,'col')
    for i=1:N_sec, sec(i).col=1; end
end
if ~isfield(sec,'shp')
    for i=1:N_sec, sec(i).shp='o'; end
end
if ~isfield(sec,'siz')
    for i=1:N_sec, sec(i).siz=zeros(N_pts,1); end
end

%% Representation of non-sector positions

sec_pos=[];
for k=1:N_sec, sec_pos=[sec_pos; sec(k).def]; end
other_pos=1:numel(v1); other_pos(sec_pos)=[];
plot(v1(other_pos),v2(other_pos),'o','MarkerSize',7,'MarkerFaceColor','w','MarkerEdgeColor','k'); 
hold on;

%% Representation of sector positions

for k=1:N_sec
    make_fig_1(v1,v2,sec(k).def,sec(k).vec,sec(k).col,sec(k).shp,sec(k).siz);
end
hold off;

%% Sub-function

function make_fig_1(v1,v2,def,vec,col,shp,siz)

bas_size=7; % basic dot size
x1=v1(def); x2=v2(def);
% if max(vec(def))>0
%     vec_n=vec(def)/max(vec(def));
% else
%     vec_n=zeros(numel(def),1);
% end
if max(vec)>0
    vec_n=vec/max(vec);
else
    vec_n=zeros(numel(def),1);
end
if abs(max(siz(def)))>0
    siz_n=3*bas_size*siz(def)/abs(max(siz(def))); % factor 3 in size at most
else
    siz_n=zeros(numel(def),1);
end

for i=1:numel(x1)
    size=max(bas_size,siz_n(i));
    plot(x1(i),x2(i),shp,'MarkerSize',size,...
         'MarkerFaceColor','w','MarkerEdgeColor','k');
    if col>=0 & col<1 & vec_n(i)>0
        [b,g,r]=hsv2rgb([col,vec_n(i),1]);
        plot(x1(i),x2(i),shp,'MarkerSize',size,...
            'MarkerFaceColor',[b g r],'MarkerEdgeColor','k');
    end
    if col<0 % grey
        plot(x1(i),x2(i),shp,'MarkerSize',size,...
            'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','k');
    end
end