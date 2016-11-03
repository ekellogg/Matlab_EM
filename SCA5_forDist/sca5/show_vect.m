function []=show_vect(V,sec,ats,option);
% usage: []=show_vect(V,sec);
%        []=show_vect(V,sec,1);
%
% Representation of a vector with colors specified in 'sec'.
% With '1' as third argument, the gradient of color are ignored.
%
% 'sec' should be a structure with fields:
% - def
% - vec
% - col
% The field 'vec' may however be omitted but the other ones are needed.

% the input ats is the position numbering list, and should be of same size
% as V.
%
% The color code is visualized with 'make_fig(0)' : the color at angle
% alpha is coded by alpha/2*PI (alpha in radian).
% For instance: 0->red, 1/3->green, 2/3-> blue.
% Admissible color codes are in [0 1] with 1 coding for white.

% Authors: Olivier Rivoire (orivoire@rockefeller.edu) and Rama Ranganathan
% 2/2010
%
%**************************************************************************

if nargin<4, option=0; end

N_pts=numel(V); N_sec=numel(sec);

% Default values for non-specified fields
if ~isfield(sec,'vec')
    for i=1:N_sec, sec(i).vec=ones(N_pts,1); end
end

%map=.8*ones(numel(sec(1).vec),3);
map=.8*ones(N_pts,3);

for s=1:N_sec
    zp=sec(s).vec/max(sec(s).vec);
    for p=1:numel(sec(s).def)
        i=sec(s).def(p);
        if sec(s).col>=0 & sec(s).col<1
            %map(i,:)=hsv2rgb([sec(s).col,zp(i),1]);
            map(i,:)=hsv2rgb([sec(s).col,zp(p),1]);
        end
        if option==1, map(i,:)=hsv2rgb([sec(s).col,1,1]); end
    end
end

%bar([V;zeros(1,numel(V))],'group');
%colormap(map);

len=numel(ats);
hold on;
for i=1:len,
    hbar(i)=bar(i,V(i));
    set(hbar(i),'FaceColor',map(i,:),'EdgeColor',[0.7 0.7 0.7]);
end;
hold off
axis([0 len+1 0 4]);grid on
set(gca,'XTick',[1:10:len]);
set(gca,'XTickLabel',ats([1:10:len]));
set(gca,'YTickLabel','');