function []=write_pmlscript(pdb_id,chain,labels_pos,sec,output)
% usages: []=write_pmlscript(pdb_id,chain,labels_pos,sec,output)
%         []=write_pmlscript(pdb_id,chain,labels_pos,sec)
%
% Writes a Pymol script to represent the sectors on a three-dimensional
% structure (pdb_id). By default the output is 'Ouputs/pmlscript_sect.pml'.
% See make_fig for the conventions concerning the colors.  Call the script
% by running pymol and executing @pmlscript.pml (with specification of the
% appropriate path).

% Authors: Olivier Rivoire (orivoire@rockefeller.edu) and Rama Ranganathan.
%
%**************************************************************************

sep='/'; if ispc, sep='\'; end % computer dependent syntax (MAC/UNIX vs PC)

N_sec=numel(sec);
if ~isfield(sec,'vec')
    for i=1:N_sec, sec(i).vec=ones(numel(sec.def),1); end
end

% Open script file
if nargin<5,
    filename=['Outputs' sep 'pmlscript_sect.pml'];
else
    filename=[output '.pml'];
end
fid=fopen(filename, 'w');

% Current directory in a format adequate for fprintf
cdir=pwd;

% Header
fprintf(fid,['# Pymol script for representing the sectors of ' pdb_id '.\n\n']);

% Show white cartoon
fprintf(fid,'delete all\n');
fprintf(fid,['load ' cdir sep 'Inputs' sep pdb_id '.pdb, main\n']);
fprintf(fid, 'hide all\n');
fprintf(fid, 'bg_color white\n'); % background color
fprintf(fid,'show cartoon');
if isempty(chain)
    fprintf(fid,'\n');
else
    fprintf(fid,[', (chain ' chain ')\n\n']);
end
fprintf(fid, 'color white\n\n');


for s=1:N_sec
    res=labels_pos(sec(s).def);
    %z=sec(s).vec(sec(s).def)/max(sec(s).vec(sec(s).def));
    z=sec(s).vec/max(sec(s).vec);
%     for i=1:numel(sec(s).def)
%         z(sec(s).def(i))=tmp(i);
%     end
    
    listsect=[];
    for i=1:numel(res)-1, listsect=[listsect res{i} ',']; end
    listsect=[listsect res{numel(res)}];
    
    % Each sector defines an object
    fprintf(fid,['create sector_' num2str(s) ', (resi ' listsect ')']); 
    if isempty(chain)
       fprintf(fid,'\n\n');
    else
        fprintf(fid,['& (chain ' chain ')\n\n']);
    end
    fprintf(fid,['show spheres, sector_' num2str(s) '\n\n']);
    
    % Each sector residue has its color
    for i=1:numel(res)
        if(z(i)>0&sec(s).col<1)
            [b,g,r]=hsv2rgb([sec(s).col,z(i),1]); 
            colname=['colres' char(res(i))];
            fprintf(fid, ['set_color ' colname ', [' num2str(b) ',' num2str(g) ',' num2str(r) ']\n']);
            fprintf(fid, ['color ' colname ', (resi ' char(res(i)) ')']);
            if isempty(chain)
                fprintf(fid,'\n\n');
            else
                fprintf(fid,['& (chain ' chain ')\n\n']);
            end
        end
    end
    fprintf(fid,['show surface, sector_' num2str(s) '\n\n']);
    fprintf(fid, 'set transparency, 0.4\n');
end

fclose(fid);
