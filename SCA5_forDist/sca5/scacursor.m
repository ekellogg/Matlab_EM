function scacursor(xresidues,yresidues,values,cluster_linewidth)
%usage: scacursor(pos_sorted,pos_sorted,C_sorted);
% Function for easy interrogation of SCA correlation matrix data, and for
% extraction of cluster composition.  Details are below.  Note that the
% position labels (xresidues and yresidues) must be cell arrays of strings.
% So typically, we pass num2cell(pos_sorted) for xresidues and yresidues.
% RR 8/2008


%============================ scacursor  v1.7 =============================
%  By: William Lane <Lucare@Lucare.com>                               
%                    Oct 2005 (Last Updated 12/27/2005)     
%==========================================================================
% Description:
%   For displaying the residue information in sca value plots.
%   For extracting the residues in a defined cluster ... creates a string
%     for easy pasting into structure programs like PyMOL.
% How to Run:
%   scacursor will attach itself to the current figure.
%   To attach during figure creation ... In the Matlab workspace type
%     imshow(scavalues);scacursor(xresidue_pos,yresidue_pos,sca_values,linewidth));
%   To attach after figure creation:
%     Select the figure you want to use by clicking it its window.
%     Then in the Matlab workspace type
%       scacursor(xresidue_pos,yresidue_pos,sca_values);
% Program Input:
%   In order to see the residue information you need to pass the 
%     residue position list to the function.
%   xresidue position = the clustered position labeles in the x direction
%   yresidue position = the clustered position labeles in the y direction
%   scavalues = the clustered sca values
%   linewidth = the width of the line for drawing boxes around clusters.
%     Default line width (if you don't specify anything) = 1.
%   To specify only one or two inputs set the others to [], except
%     linewidth which you should leave blank if you do not wish to specify.
%   This program can be called with no inputs if you wish ... though there
%     is really no point to that.
% Example Use:
%   [p_pos,l_pos,sort_pos,sorted]=rr_cluster_2(cmr_beta,pos_beta,3.5,'jet',0);
%   scacursor(pos_beta(sort_pos),pos_beta(sort_pos),sorted);
% Program Use:
%   With the figure window active type 'h' to toggle the help box on|off.
%   With the figure window active type 'd' to toggle on|off the disaplay 
%     that follows the cursor.
%   With the figure window active type 'c' to draw a rectangle around the
%     cluster you wish to get residue info from ... for cutting and pasting
%     into a structure program like PyMOL. The cluster will be defined by
%     those pixels within the box and those that the box's line touch.
%   With the figure window active type 'q' to exit scacursor.
%   Note: you can use the zoom and pan tools. Just unselect the active tool 
%     to get back to the scacursor.
%==========================================================================
%     Version History
%-----------------------------------------------
%     v0.1 - basic cursor movment working with clicking for zoom and mouse
%            panning.
%     v0.2 - removed clicking for zoom and panning (it is easier to use the
%            ones in the figure menu)
%     v0.3 - figured out how to get the integer pixel number from the 
%            cursor position
%     v0.4 - added function inputs for determining residue and sca values
%     v0.5 - added cursor movement with displayed 
%            values following movment
%     v0.6 - added keyboard shorcut for toggling cursor following display
%     v0.7 - basic rubber band functionality
%     v0.8 - fixed residue display ... used {} instead of () since the 
%            varible was of type cell.
%     v0.9 - fixed cluster rectangle drawing on windows xp
%     v1.0 - added displaying of cluster info in matlab workspace ... 
%            for cutting and pasting to structure program
%     v1.1 - fixed reseting key and mouse movement binding for figure 
%            on scacusor exit
%     v1.2 - added sca value display
%     v1.3 - added help box
%     v1.4 - added user specified line width for cluster drawing box
%     v1.5 - changed residue display back as in v0.8 used {} instead of 
%            () since the varible was of type cell.
%     v1.6 - fixed value display ... x and y were backwards (ie
%            value = M(y,x) and not value = M(x,y).
%     v1.7 - added left and right arrow for changing the Max image value.
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
% 
% scacursor
% main/primary function 
% get figure handles, create the gui elements, etc ...
%--------------------------------------------------------------------------

    % set to show debug info. (0|1) - 0: do not show, 1: show
    scacursorparams.debug = 0;

    % get the function inputs and create empty varible it not passed   
    if nargin < 1, xresidues = []; end
    if nargin < 2, yresidues = []; end
    if nargin < 3, values = []; end
    if nargin < 4, cluster_linewidth = 1; end
    scacursorparams.xresidues = xresidues;
    scacursorparams.yresidues = yresidues;
    scacursorparams.values = values;
    scacursorparams.cluster_linewidth = cluster_linewidth;
    
    % make sure that there is a current figure
    currfig = findobj('type','figure');
    if isempty(currfig)
        beep;
        h=warndlg('There are no figures open!');
        uiwait(h);
        return
    end

    % make sure that there is a an axis on this figure
    currfig = currfig(1);
    figure(currfig);
    scacursorparams.currax = findobj(currfig,'type','axes');
    scacursorparams.currax = get(currfig,'currentaxes');
    if isempty(scacursorparams.currax)
        error('The current figure contains no axes!');
    end

    % make sure that there is an image in this axis
    scacursorparams.currobj = findobj(currfig,'type','image');
    if isempty(scacursorparams.currobj)
        error('The current axis doesn''t appear to contain any valid images.');
    end
    
    % make the pointer a crosshair
    scacursorparams.CLim_max = 2;
    set(scacursorparams.currax, 'CLim',[0,scacursorparams.CLim_max]);
    
    % capture mouse movement ... will call positionfcn everytime the mouse moves
    %setappdata(currfig, 'CLim', 0);    
    
    % make the pointer a crosshair
    set(currfig, 'Pointer', 'crosshair','units','normalized');
    
    % get the current axis
    currax = get(scacursorparams.currax);

    % capture key presses ... will call kpfcn everytime a key is pressed
    set(currfig, 'KeyPressFcn',@kpfcn);   
    
    % capture mouse movement ... will call positionfcn everytime the mouse moves
    setappdata(currfig, 'positionfcnhandle', @positionfcn);
    set(gcf, 'WindowButtonMotionFcn','feval(getappdata(gcf,''positionfcnhandle''));');

    % create the help text box
    scacursorparams.help_visible = 'off';
    scacursorparams.dispbox2 = text('Parent', scacursorparams.currax, 'HorizontalAlignment', 'left', ...
        'Visible', scacursorparams.help_visible, ...
        'Color', 'y', 'Tag', 'HELP',...
        'backgroundcolor','k', 'Position', [0 100],...
        'FontSize', get(scacursorparams.currax, 'FontSize') / 1);  
    set(scacursorparams.dispbox2,'string',sprintf('Type ''h'' to hide and show this text box.\nHit ''left arrow key'' to decrease the color range.\nHit ''left arrow key'' to increase the color range.Type ''d'' to toggle the cursor display on|off.\nType ''c'' then left click start point and then left click end point to select a cluster to display residue info about.\n    The cluster will be defined by those pixels within the box and those that the box''s line touch.\nType ''q'' to exit scacursor\nNote: you can use the zoom and pan tools. Just unselect them to get back to the scacursor.'));
    
    % create the text box that moves with the cursor and displays the values
    % key is to make its Parent the current axis ... so that its coords are relative to the image
    scacursorparams.show_cursor_values_visible = 'on';
    scacursorparams.dispbox3 = text('Parent', scacursorparams.currax, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'Visible', scacursorparams.show_cursor_values_visible, ...
        'Color', 'y', 'Tag', 'CURSOR_VALUES',...
        'backgroundcolor','k',...
        'FontSize', get(scacursorparams.currax, 'FontSize') / 1);
    
    %scacursorparams.dispbox4 = text('Parent', scacursorparams.currax, 'HorizontalAlignment', 'left', ...
    %    'Visible', 'on', ...
    %    'String','CLim', 'Color', 'y', 'Tag', 'HELP',...
    %    'backgroundcolor','k', 'Position', [0 0],...
    %    'FontSize', get(scacursorparams.currax, 'FontSize') / 1.2);  
    %set(scacursorparams.dispbox4,'string',sprintf('%s',num2str(scacursorparams.CLim_max)));
    %scacursorparams.slider = uicontrol('Style','slider','Max', 4, 'Min', 0, 'callback',...
    %    ['scacursorparams = getappdata(gcf,''scacursorparams'');set(scacursorparams.currax, ''CLim'',[0,get(scacursorparams.slider,''Value'')]);fprintf(''%s\n'',num2str(get(scacursorparams.slider,''Value'')));']);
    
    % turn on the figure toolbar ... so that we can have quick access to the zoom controls
    set(gcf, 'Toolbar','figure');
    
    % set rubber band box as not being displayed ... this way we can turn it on later and have positionfcn only deals with it when it is being used.
    scacursorparams.displaying_rubberband_box = 0;
    
    % save the current structure variable scacursorparams to this figure.
    % this allows the functions and the callbacks to talk to each other
    setappdata(gcf,'scacursorparams',scacursorparams);
%==========================================================================

%==========================================================================
% kpfcn
%   called upon a key press ...
%   takes care of keypress:
%    - 'h' then toggle the help box on|off
%    - 'd' then toggle on|off the value display that follows the cursor
%    - 'c' rubberband around the cluster and writeou PyMOL code
%    - 'q' quit the program
%--------------------------------------------------------------------------
function kpfcn(src,evnt)
    % read in the saved structure scacursorparams    
    scacursorparams = getappdata(gcf,'scacursorparams');   
    
    % debug info
    if scacursorparams.debug == 1, fprintf('---------- In key_press - %s\n', evnt.Character); end
  
    % evaluate the the keypress ...
    
    % if keypress is 'h' then toggle on|off the help box
    if evnt.Character == 'h'
        if   strcmp(scacursorparams.help_visible, 'on')
             scacursorparams.help_visible = 'off';
        else
            scacursorparams.help_visible = 'on';
        end
        set(findobj('Tag','HELP'),'Visible',scacursorparams.help_visible);
    end

    % change the color range of the image using the left and right arrow keys
    if evnt.Character == 28     % left arrow key
       if scacursorparams.CLim_max > 0.1    % max can't be lower than our min
           scacursorparams.CLim_max = scacursorparams.CLim_max - 0.1;
           set(scacursorparams.currax, 'CLim',[0,scacursorparams.CLim_max]);
           %set(scacursorparams.dispbox4,'string',sprintf('%s',num2str(scacursorparams.CLim_max)));
           set(scacursorparams.dispbox3,'string',sprintf('  image range=[0 %s]',num2str(scacursorparams.CLim_max)));
           % disp(scacursorparams.CLim_max);
       end
    end
    if evnt.Character == 29     % right arrow key
       scacursorparams.CLim_max = scacursorparams.CLim_max + 0.1;
       set(scacursorparams.currax, 'CLim',[0,scacursorparams.CLim_max]);
       set(scacursorparams.dispbox3,'string',sprintf('  image range=[0, %s]',num2str(scacursorparams.CLim_max)));
       %set(scacursorparams.dispbox4,'string',sprintf('%s',num2str(scacursorparams.CLim_max)));
       %disp(scacursorparams.CLim_max);
    end
    
    % if keypress is 'd' then toggle on|off the value display that follows the cursor
    if evnt.Character == 'd'
        if   strcmp(scacursorparams.show_cursor_values_visible, 'on')
             scacursorparams.show_cursor_values_visible = 'off';
        else
            scacursorparams.show_cursor_values_visible = 'on';
        end

        % only turn the cursor values text back on right away if we are in the image
        % get the x, y  positions
        posn = get(scacursorparams.currax,'currentpoint');
        posn = posn(1,:);
        x = posn(1,1);
        y = posn(1,2);   
        x_int = ceil(x-0.5);        % this is a integer corresponding to the box pixel number in the x
        y_int = ceil(y-0.5);        % this is a integer corresponding to the box pixel number in the y
        % get the current image boundary values
        xlim_current = xlim(scacursorparams.currax);
        ylim_current = ylim(scacursorparams.currax);
        % check if the cursor is within the current image and display position
        if x_int >= xlim_current(1) & x_int <= xlim_current(2) & x_int > 0 & y_int >= ylim_current(1) & y_int <= ylim_current(2) & y_int > 0 
            set(findobj('Tag','CURSOR_VALUES'),'Visible',scacursorparams.show_cursor_values_visible);
        end
    end

    % if keypress is 'c' then draw a rubberband around the cluster and writeout the text to create object in PyMOL
    if evnt.Character == 'c'
        % debug info
        if scacursorparams.debug == 1, sprintf('---------- In rubberband\n'); end
        
        % Get current user data ... using UserData to share info instead of setappdata ... to keep the two seperate.
        cudata=get(gcf,'UserData'); 
        hold on;
        
        % Wait for for mouse click ... and get start point
        k=waitforbuttonpress;
        p1=get(scacursorparams.currax,'CurrentPoint');       %get starting point
        p1=p1(1,1:2);                   % set p1 to x,y
        lh=plot(p1(1),p1(2),'LineWidth', scacursorparams.cluster_linewidth,'Color','g');      %plot starting point
        % debug info
        if scacursorparams.debug == 1, sprintf('p1(1)=%3.2f     p1(2)=%3.2f\n',p1(1),p1(2)); end
        udata.p1=p1;
        udata.lh=lh;
        set(gcf,'UserData',udata);
        scacursorparams.displaying_rubberband_box = 1;
        setappdata(gcf,'scacursorparams',scacursorparams);
        
        % Wait for for mouse click ... and get end point
        k=waitforbuttonpress;
        p2=get(scacursorparams.currax,'Currentpoint');       %get end point
        p2=p2(1,1:2);                   % set p1 to x,y

        set(gcf,'UserData',cudata); %reset UserData, etc..
        delete(lh);
        % debug info
        if scacursorparams.debug == 1, sprintf('p2(1)=%3.2f     p2(2)=%3.2f\n',p2(1),p2(2)); end
        
        scacursorparams.displaying_rubberband_box = 0;
        setappdata(gcf,'scacursorparams',scacursorparams);

        p11_int = ceil(p1(1)-0.5);        % this is a integer corresponding to the box pixel number in the x
        p12_int = ceil(p1(2)-0.5);        % this is a integer corresponding to the box pixel number in the y
        p21_int = ceil(p2(1)-0.5);        % this is a integer corresponding to the box pixel number in the x
        p22_int = ceil(p2(2)-0.5);        % this is a integer corresponding to the box pixel number in the y
        if scacursorparams.debug == 1, sprintf('p11_int=%d    p12_int=%d\n',p11_int,p12_int); end
        if scacursorparams.debug == 1, sprintf('p21_int=%d    p22_int=%d\n',p21_int,p22_int); end        
        
        % display the cluster residues ... 
        % cut and pasted into a structure program for easy display of clustered residues.
        % get the current image boundary values
        xlim_current = xlim(scacursorparams.currax);
        ylim_current = ylim(scacursorparams.currax);
        % loop through the cluster in the x dir
        p_start = min(p11_int,p21_int);
        p_end = max(p11_int,p21_int);
        cluster_string_to_print = '';
        for p_index = p_start:1:p_end
            %disp(p_index);
            % check if the cursor is within the current image and display position
            if p_index >= xlim_current(1) & p_index <= xlim_current(2) & p_index > 0
                if ~(isempty(scacursorparams.xresidues)) & p_index <= size(scacursorparams.xresidues,2)
                    %disp(scacursorparams.xresidues{p_index});       % this is a cell array ... using the curly bracket gets rid of the ''
                    if strcmp(cluster_string_to_print,'')
                        cluster_string_to_print = strcat(cluster_string_to_print,num2str(scacursorparams.xresidues{p_index}));
                    else
                        cluster_string_to_print = strcat(cluster_string_to_print,'+',num2str(scacursorparams.xresidues{p_index}));
                    end
                end
            end
        end
        disp('---------------------------------');
        if ~(strcmp(cluster_string_to_print,''))
            disp('residues in the cluster''s x dir ... for pasting into structure program');
            disp(cluster_string_to_print);
        else
            disp('no residues in the x dir ... did you input the position array');
        end
        disp('---------------------------------');
        % loop through the cluster in the y dir
        p_start = min(p12_int,p22_int);
        p_end = max(p12_int,p22_int);
        cluster_string_to_print = '';
        for p_index = p_start:1:p_end
            %disp(p_index);
            % check if the cursor is within the current image and display position
            if p_index >= ylim_current(1) & p_index <= ylim_current(2) & p_index > 0
                if ~(isempty(scacursorparams.yresidues)) & p_index <= size(scacursorparams.yresidues,2)
                    if strcmp(cluster_string_to_print,'')
                        cluster_string_to_print = strcat(cluster_string_to_print,num2str(scacursorparams.yresidues{p_index}));  % we need curly brackets since this is a cell variable
                    else
                        cluster_string_to_print = strcat(cluster_string_to_print,'+',num2str(scacursorparams.yresidues{p_index}));   % we need curly brackets since this is a cell variable
                    end
                end
            end
        end
        if ~(strcmp(cluster_string_to_print,''))
            disp('residues in the cluster''s y dir ... for pasting into structure program');
            disp(cluster_string_to_print);
        else
            disp('no residues in the y dir ... did you input the position array');
        end
            disp('---------------------------------');
    end
    
    % if keypress is 'q' then quit
    if evnt.Character == 'q'
        set(gcf,'KeyPressFcn','');      % stop evaluting kep presses
        set(gcf, 'WindowButtonMotionFcn','');   % stop evaluating mouse movement
        delete(scacursorparams.dispbox2);       % delete help box
        delete(scacursorparams.dispbox3);       % delete cursor value box
        %delete(scacursorparams.dispbox4);       % delete color range (CLim) value box
        %delete(scacursorparams.slider);         % delete the slider
    end
    
    % save the current structure variable scacursorparams to this figure.
    % this allows the functions and the callbacks to talk to each other
    setappdata(gcf,'scacursorparams',scacursorparams);
return
%==========================================================================

%==========================================================================
% positionfcn
%   called upon mouse movement
%   takes care of changing things based on mouse movment:
%   (1) getting the cursor position and updating display values
%   (2) draw the rubber band if applicable
%--------------------------------------------------------------------------
function positionfcn
    % read in the saved structure scacursorparams
    scacursorparams = getappdata(gcf,'scacursorparams');
    
    % debug info
    if scacursorparams.debug == 1, sprintf('---------- In positionfcn\n'); end

    % get the x, y  positions
    posn = get(scacursorparams.currax,'currentpoint');
    posn = posn(1,:);
    x = posn(1,1);
    y = posn(1,2);   
    x_int = ceil(x-0.5);        % this is a integer corresponding to the box pixel number in the x
    y_int = ceil(y-0.5);        % this is a integer corresponding to the box pixel number in the y
        
    % debug info
    if scacursorparams.debug == 1, sprintf('x=%3.0f, y=%3.0f\n',x,y); end

    % get the current image boundary values
    xlim_current = xlim(scacursorparams.currax);
    ylim_current = ylim(scacursorparams.currax);
    
    % debug info
    if scacursorparams.debug == 1, sprintf('xlim=%3.0f, ylim=%3.0f\n',xlim_current(1),ylim_current(1)); fprintf('xlim=%3.0f, ylim=%3.0f\n',xlim_current(2),ylim_current(2)); end
    
    % check if the cursor is within the current image and display position
    if x_int >= xlim_current(1) & x_int <= xlim_current(2) & x_int > 0 & y_int >= ylim_current(1) & y_int <= ylim_current(2) & y_int > 0 
        if ~(isempty(scacursorparams.xresidues)) & x_int <= size(scacursorparams.xresidues,2)
            xresidue_to_print = scacursorparams.xresidues{x_int};      % we need curly brackets since this is a cell variable
        else
            xresidue_to_print = 'n/a';
        end
        if ~(isempty(scacursorparams.yresidues)) & y_int <= size(scacursorparams.yresidues,2)
            yresidue_to_print = scacursorparams.yresidues{y_int};       % we need curly brackets since this is a cell variable
        else
            yresidue_to_print = 'n/a';
        end
        if ~(isempty(scacursorparams.values)) & x_int <= size(scacursorparams.values,2) & y_int <= size(scacursorparams.values,1)
            sca_value_to_print = sprintf('%3.2f',scacursorparams.values(y_int,x_int));
        else
            sca_value_to_print = 'n/a';
        end
        % update residue display text
        set(scacursorparams.dispbox3,'position', [x y] ,'String',sprintf('x = %s        y = %s\n  value = %s    [%s]',num2str(xresidue_to_print),num2str(yresidue_to_print),num2str(sca_value_to_print),num2str(scacursorparams.CLim_max)));
        if strcmp(scacursorparams.show_cursor_values_visible, 'on')
                set(findobj('Tag','CURSOR_VALUES'),'Visible','on');
        end
    else
        if strcmp(scacursorparams.show_cursor_values_visible, 'on')
                set(findobj('Tag','CURSOR_VALUES'),'Visible','off');
        end
    end

    if scacursorparams.displaying_rubberband_box == 1
        utemp=get(gcf,'UserData');
        ptemp=get(scacursorparams.currax,'CurrentPoint');
        ptemp=ptemp(1,1:2);
        % draw box
        set(utemp.lh,'XData',[ptemp(1),ptemp(1),utemp.p1(1),utemp.p1(1),ptemp(1)],'YData',[ptemp(2),utemp.p1(2),utemp.p1(2),ptemp(2),ptemp(2)]);
        drawnow;    % required to show the rectangle on windows xp
    end
    
    % debug info
    if scacursorparams.debug == 1, sprintf('----------\n'); end
return
%==========================================================================



