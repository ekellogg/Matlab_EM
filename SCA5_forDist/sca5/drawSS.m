function drawSS(pdb,ats,ss,sec,colors)
% usage: drawSS(pdb,alignment_to_structure)
%        drawSS(pdb, alignment_to_structure,ss)
%        drawSS(pdb, alignment_to_structure,ss,sec)
%        drawSS(pdb, alignment_to_structure,ss,sec,['b' 'r' 'g'])

% This function takes in a pdb file and a list of positions (ats) and
% outputs a graph showing the secondary structure pattern of the protein
% trimmed to the list of selected positions. Call the function after
% creating a figure.  For example:
%
% figure(100);drawSS(pdb,ats,ss,sec);
%
% The optional input "ss" is a structure with fields ss.Helix and ss.Sheet,
% and is specified if the secondary structure to be shown is a subset of
% the annotations in the pdb file.  ss should contain the number of each
% secondary structure element to be represented.  The actual starting and
% ending positions are determined from the pdb file.  For example:
%
% ss = 
%    Sheet: [1 2 3 6 7] 
%    Helix: [2 3]
%
% indicates and Sheets 1,2,3,6, and 7 should be shown and helices 2 and 3
% should be shown. The optional input variable "sec" is a structure that
% includes a field called sec.def, which contains the indices for positions
% in ats that comprise each sector. The first is taken as the blue sector,
% and if they exist, the second as the red sector, the third as the green
% sector, and fourth as yellow, and the fifth as cyan.  IF a different
% coloring scheme is desired, pass this in the string array "colors".  For
% example, if the order should be blue, red, and green"
%
% colors = ['b' 'r' 'g']
%
% The position list "ats" can be a numeric vector, a numeric cell array, or
% a cell array of strings.  The latter is useful for some pdbs that contain
% position values that include characters (e.g. position 221A).
%
% Author: Rama Ranganathan
% 12/2009
% copyright Rama Ranganathan, 2009-2010

% preliminiaries
sec_index = 1;

if nargin<4
    sec_index = 0;
end
if nargin<5
    colors=['b' 'r' 'g' 'y' 'c'];
end

% convert all ats formats to a cell array of strings

if ~iscellstr(ats)
    if iscell(ats)
        ats=cell2mat(ats);
    end
    if isnumeric(ats)
        if size(ats,2)>1;ats=ats';end;
        ats=cellstr(num2str(ats));
    end
end
len=numel(ats); 


% draw an initial line representing the whole protein
%clf
h1=plot([1:len],repmat(0,1,len),'k','LineWidth',5);grid on
axis([0 len+1 -0.5 0.5]);
set(gca,'XTick',[1:10:len]);
set(gca,'XTickLabel',ats([1:10:len]));
set(gca,'YTickLabel','');

% draw the Sheets
jnk=squeeze(struct2cell(pdb.Sheet));jnk=jnk([6 10],:);  
if nargin>=3; jnk=jnk(:,ss.Sheet); end;
for i=1:size(jnk,2)
    if ~isempty(strmatch(num2str(jnk{1,i}),ats,'exact'))
        sheet_start=strmatch(num2str(jnk{1,i}),ats,'exact');
        sheet_end=strmatch(num2str(jnk{2,i}),ats,'exact');
        hold on;plot([sheet_start:sheet_end-2],repmat(0,1,numel([sheet_start:sheet_end-2])),'k','LineWidth',10);
        hold on;plot(sheet_end-1,0,'>k','MarkerSize',15,'MarkerFaceColor','k');
    end
end

%draw the helices
jnk=squeeze(struct2cell(pdb.Helix));jnk=jnk([5 9],:);  
if nargin>=3; jnk=jnk(:,ss.Helix); end;
for i=1:size(jnk,2)
    if ~isempty(strmatch(num2str(jnk{1,i}),ats,'exact'))
        helix_start=strmatch(num2str(jnk{1,i}),ats,'exact');
        helix_end=strmatch(num2str(jnk{2,i}),ats,'exact');
        hold on;plot([helix_start:helix_end],repmat(0,1,numel([helix_start:helix_end])),'k','LineWidth',10);
    end
end

% add the sectors as colored squares if sec is specified
if sec_index~=0
    numsec=numel(sec);
    hold on
    for i=1:numsec
        scatter(sec(i).def,repmat(0,1,numel(sec(i).def)),['s' colors(i)],'SizeData',25,'MarkerFaceColor',colors(i));
    end
    hold off
end

end

