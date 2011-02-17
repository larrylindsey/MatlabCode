function varargout = gui(varargin)
% Graphical User Interface for analyzing automatic STEM segmentations.
% Message from Matlab:
%      "gui", by itself, creates a new GUI or raises the existing
%        singleton.
%      h = gui returns the handle to a new GUI or the handle to
%        the existing singleton.
%      gui('CALLBACK',hObject,eventData,handles,...) calls the local
%        function named CALLBACK in gui.m with the given input arguments.
%      gui('Property','Value',...) creates a new GUI or raises the
%        existing singleton*.  Starting from the left, property value 
%        pairs are applied to the GUI before gui_OpeningFunction gets 
%        called.  An unrecognized property name or invalid value makes 
%        property application stop.  All inputs are passed to 
%        gui_OpeningFcn via varargin.
%
% GLOBAL STARTUP FLAGS:
%  flNoBuffLoad     disables default loading buffer file at startup
%  flMonitor        disables backup/buffer file updates, used for
%                  inspecting someones work in progress
%
% v7.5  written by Y.Mishchenko  JUL 2006   Chklovskii Lab   CSHL


% updates LOG
% started from v. 7.1 OCT 15, 2006
% v7.50 JUN 10
% * fixed bwlabel-once to repeat bwlabel after draw-boxes had been added
% * made gui to work with file-by-file format of data
% * had to change drawing editions are stored in data.tmp{20} to deal with
%   file-by-file data format
% v7.40 MAY 21
% * added reduced image routine to speed up screen updates for large images
% * made bwlabel current section only once when quick-draw for faster
%   drawing, this may be buggy needs attention in the future
% * made selection input area loose focus on enter
% v7.33 MAR 15
% * added 3D thumbnail function
% * changes/improvements for in-memory-buffer
% v7.32 MAR 05
% * added reminder before exit
% * changed default hotkey assignment
% * changed gui.backup.mat to gui.backup###.mat
% * changed save-method for buff for -v6 == no compressing, faster
% v7.31 FEB 13
% * fix bug with proof from gui.backup.mat been overwritten 
%   from gui.dataXXX.mat.
% v7.3 JAN 29
% * fix the bug with push id, leading to data.cmap causing out-of-bounds
%   index after pusid
% * modified the way log of drawing alterations is updated, to save update
%   time when there is large amount of drawings
% * changed to imcomplement in special overlay mode
% * introduced memory buffer to remove time overhaul during frequent
%   switching between adjacent sections, to facilitate binding corrections
%   use global flUseMemBuffer=1 for that; be careful with out-of-memory's
% * flMonitor=1 now does not update buffer during crossovers, prints
%   warnings, flMonitor=1 now still allows to manually save backup
% * proof notes now may be saved in "wrong" substack if they were
%   initialized in that substack
% * proof list now will NOT be updated in the "wrong substack" upon call
%   save-proof-note function
% * when performing binding subseries, set global flMonitor to 1, 
%   this will remove delays due to buffer file update when 
%   alternating between sections, then remember to manually save updates 
%   for target section when binding correction is complete
% v.7.2 DEC 28
% * fixed bug with out-of-range colormaping when crossing to adjacent
%   subseries
% * changed shading from inverse to direct percentage-brighter (>1)
% * changed imagesc to imshow and made changes to graphic dataset
%   generation cutting about 30% time from color screen updates
% * experiment alternative color mode
% * experiment shadow mode, may supply autosegmentation reduced to 1pxl
%   separators in debug.ind2 for use with shadow mode
%   code to produce reduced segmentation is this (assuming cat is uint16)
%   ind2={}; 
%   for i=1:length(cat) 
%   X=watershed(imimposemin(imcomplement(al{i}),cat{i}>0))==0;
%   Y=zeros(size(X),'uint16'); Y(:)=65535; Y(X)=0; 
%   ind2{i}=imreconstruct(uint16(cat{i}),Y)); 
%   fprintf('.'); end; fprintf('\n');
% * fixed hotkey-handler missing on "play" button
% * fixed selection highlight so it works faster for smaller FOVs
% v.7.1 OCT 15
% * added dialog boxes when asking to load backup & for error message
% * added global flag to disable loading of buff file flNoBuffLoad
% * added global flag for 'monitor' mode - no saves for buffs flMonitor
% * added proof.ttime array to track timing of the corrections
% * added verification al is uint8
% * changed double-move-to-cross to triple-move-to-cross
% * changed focus shift to reset - fixes hotkey disabling at top section
% * fixed problem with resizing when switching to major mode
% * added sorting options to proof-list
% * altered next/prev button to advance through slice-specific selection in
%   proof-list, correspondingly, removed references to "super" mode as it
%   is now implementable via "By Slice" filter in proof-list
% * added size-reduction for cat depending on largest object ID within
%   startgui to reduce memory requirements
% * changes made to loading backups allowing loading backup for not active 
%   miniseries via requesting "load backup" while in that miniseries


% FROM V7 DISCOUNTINUED BACKWARD COMPATIBILITY FOR EARLIER
% SEGMENTED STACKS HARRIS S01, S06; KNOTT S02, S08, USE V5 FOR THOSE


% THIS IS OUTDATED:
% -----------------------------------------------------------------------
% Input stack & debug data & stuff are as follows:
% stack of original grayscale images - 'al'
% segmentation in label-array format - 'cat'
% segmentation info  - 'debug' structure; includes
%
% initial information:
%  .slc     - original 2D labels
%  .cat0    - forward pass ['working'] labels
%  .mskE    - logical array of 'exposed' mask == edges
% processed information
%  .logs    - structarray of log records for events
%  .map     - correspondence cat0 -> cat map
%  .ide     - [debug.logs.event] for gui
%  .idd     - [debug.logs.outgoing] for gui
%  .idk     - [debug.logs.slice]  for gui
%  .idt     - map transformed [debug.logs.tag]
%  .esc     - list of escaping (major) clusters for gui
%  .sstats  - array with statistics for each cluster for gui, 
%     these are max amongst all slices containing the fragment:
%      [gray mean [0-1], anchors mean, majoraxis, minoraxis, area/1000,... 
%             extent in terms of slices it goes through]
% -----------------------------------------------------------------------

% THIS IS OUTDATED:
% -----------------------------------------------------------------------
% output data are as follows:
% proof     - proof reading structure including
%  .pmap    - final cat0 -> cat labels assignments
%  .tmap    - process type ids:
%               1 - none; 2 - dendrite; 3 - axon; 4 - glia; 
%               5 - extracell space; 8 - detached spine; 11 - other
%  .notes   - structure containing operator notes; each contains
%   .logid  - id of the log record being examined
%   .concl  - operator conclusion ==1*corrected+2*notsure+4*other
%   .stats  - additional stats info about corrections made by operator
%   .note   - text note made by operator
%  .v       - proof version
%  .ttime   - proof-timing information
% other fields


% ======================================================================
% Edit the above text to modify the response to help gui
% Last Modified by GUIDE v2.5 14-Mar-2007 23:34:56
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% ======================================================================


% ======================================================================
%                           INITIALIZATION  
% ======================================================================
% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% global input data
global al cat debug proof

% internal shared data - primary storage
global data

% shared tmp variables
global img idx tmp tmp3

% controller for noBuffLoad
global flNoBuffLoad

% controller for monitoring
global flMonitor

% check if first use, used in reinitialization in utdmove & savebackup
firstuse=isempty(eventdata);

if(firstuse)
    fprintf('Version 7.50\n\n');
    fprintf('NOTE: gui is using GLOBAL variables to share information,\n');
    fprintf('please do not ''clear all'' when gui is running...\n\n');
    fprintf('gui is initializing its data environments...\n');
    
    global buffer
    buffer={};    
end


% -----------------------------------------------------------------------
% CHECK CRITICAL PARAMETERS:
% make sure that we have our slices
if(isempty(al))
    flquit=0;
    if(~isfield(debug,'prefix') | isempty(debug.prefix)) flquit=1;
    else    % initialize loading from HD
        flist=dir([debug.prefix,'*','.mat']);
        if(isempty(flist)) flquit=1;
        else
            kk=[];
            for k=1:length(flist)
                sname=flist(k).name;
                sname=sname(end-4:-1:1);
                sname=sname(1:find(sname=='.',1));
                sname=sname(end-1:-1:1);
                kk=[kk,str2num(sname)];
            end
            
            if(isempty(kk)) flquit=1;
            else
                al=cell(1,max(kk));
                cat=cell(1,max(kk));
                
                sname=sprintf('%s.001.mat',debug.prefix);
                
                if(~exist(sname)) flquit=1; 
                else
                    wc=load(sname,'al','cat'); al{1}=wc.al; cat{1}=wc.cat;
                end
            end
        end
    end                

    if(flquit)
        errordlg('ERROR WHILE LOADING DATA!!!'); 
        error('ERROR WHILE LOADING DATA'); 
    end
end

% make sure al is in cell formats
if(~iscell(al)) 
    fprintf('GUI v7 uses EM data in cell format, attempting to convert...\n');
    al1=cell(1,size(al,3));
    for k=1:size(al,3) al1{k}=al(:,:,k); end
    al=al1; clear al1;
    fprintf('convert EM stack successful\n');
end

% make sure our slices are uint8
for k=1:length(al) al{k}=im2uint8(al{k}); end

% if running as viewer for al
if(isempty(cat))
    cat=cell(size(al));
    for k=1:length(al)
        cat{k}=false(size(al{k})); cat{k}(1,1)=1;
    end
end

% make sure cat is in cell format
if(~iscell(cat)) 
    fprintf('GUI v7 uses segmentation in cell format, attempting to convert...\n');
    cat1=cell(1,size(cat,3));
    for k=1:size(cat,3) cat1{k}=cat(:,:,k); end
    cat=cat1; clear cat1;
    fprintf('convert segmentation successful\n');
end



% INITIALIZe VARIABLES

% *************** JUMPERS *******************
% set this flag to 1 to force GUI to build anchors and membranes
% whenever absent from supplied 'debug' var
force_fields=0;
% set this flag to 1 to allow 'gather' mode in GUI
allowgather=1;
% use ordering of proofing-list
prooforder=1;

% default threshold for anchors
bgthr=0.7;
% size of preallocated proof-notes array
maxnote=10^4;

if(firstuse)
data.version=7;                     % version of the proof file 
                                    % set 3 for s01.KHarris                                        
data.allowgather=allowgather;       % allow gather function?
data.edgemax=6;                     % highest gradation in edges map 


% PREFERENCES
data.pref.w=100;                    % how many pxl additionally allow 
                                    % around a fragment for closeup
data.pref.shade=[1.5,1.25];        % how much insignificant segments are 
                                    % shaded in gui: [primary,alternative]
data.pref.clmp='jet';               % colormap to use when drawing
data.pref.edrange=2;                % pen size for editing membranes
data.pref.light=1;                  % use light mode?
data.pref.selmode=0;                % selection mode?
data.pref.dbg=0;                    % show debug records?
data.pref.objdskt=[2,4,11];         % object types not discounted thr-mode
% these are preset hotkeys order, do not alter in mainKeyDownFcn
data.pref.keydscrp={'up';'down';'left';'right';'zoom in';...
    'zoom out';'move up';'move down';'flip';'cycle';'overlay';...
    'next';'previous';'quick draw';'quick link';'redraw';'undo';...
    'quick del';'delete';'draw';'recenter'};
% hotkeys in the order set above
data.pref.keys=[30,31,28,29,61,45,101,119,113,99,111,110,112,115,...
    114,108,26,97,102,100,120];
data.pref.watershed=1;              % whether watershed is for single section?
data.pref.undocleanup=0;            % whether cleanup undo after each obj?
data.pref.selmajor=1;               % select only major on box-select?
data.pref.altsel={'normal','extend','alt'};
                                    % click-selection preferences
data.pref.todolength=50;            % size of page in to-do list

data.pref.dftcor=2;                 % default corrected state - here prev
data.pref.dftcklink=1;              % default quick-link state
data.pref.dftpen=2;                 % default drawing pen

data.pref.con=4;                    % connectivity used in quick-edit modes
data.pref.delbuf=0;                 % delete cross-stack buffers on exit?
data.pref.ccdouble=1;               % press move triple before changing stack?

% update preferences from file if any
if(exist('pref.mat')) load pref.mat; data.pref=prefs;  end
end

% identify size of stack and position of first section
% ment to be used with "sliding" substack, but no implementation
data.csize=Inf;
data.cpos=1;


% other internal parameters
data.shape=size(cat{1});                % image size
data.smax=min(length(al),length(cat));  % extent in 3rd dim
data.v1zoom=[];                         % processed image window
data.v2zoom=[];                         % processed selection window

% max object # in cat
if(isfield(debug,'unique')) data.mmax=double(max(debug.unique));
else data.mmax=0; end
for i=1:length(cat) 
    if(~isempty(cat{i}))
        data.mmax=max(data.mmax,max(cat{i}(:))); 
    end
end
data.mmax=double(max(1,data.mmax));

% data-size reduction for cat
if(data.mmax==1) 
    for i=1:length(cat) cat{i}=logical(cat{i}); end
elseif(data.mmax<2^8)
    for i=1:length(cat) cat{i}=uint8(cat{i}); end
elseif(data.mmax<2^16)
    for i=1:length(cat) cat{i}=uint16(cat{i}); end
end

% tracking variables
data.ids=[];            % global ids of selected objects
data.logid=0;           % id of log record selected
data.nclst=0;           % id of last selected process

if(firstuse)
% current zoom on the stack
data.zoom=[1,size(al{1},1),1,size(al{1},2)];
data.tmp=cell(20,1);    % temporary array, currently as follows:
                        % .tmp{1} used to store first overlay member
                        % .tmp{2} logids selection for '...'-mode
                        % .tmp{3} shift-click buffer
                        % .tmp{4} holds first segment for link
                        % .tmp{5} holds slice # in proofing process
                        % .tmp{6} holds second segment for link
                        % .tmp{7} vacated
                        % .tmp{8} vacated
                        % .tmp{9} stores items in prooflist that corrected
                        % .tmp{10} vacated
                        % .tmp{11} ids of all 'notsure' recs in
                        %          'supervisor' reading mode
                        % .tmp{12} vacated
                        % .tmp{13} vacated
                        % .tmp{14} vacated
                        % .tmp{15} vacated
                        % .tmp{16} vacated
                        % .tmp{17} holds bwlabel for quick-draw edits
                        % .tmp{18} holds current process id
                        % .tmp{19} holds list of all unique objects in cat
                        % .tmp{20} list of drawing alterations to cat

% define selection in dropboxes
if(~data.pref.light)
    data.dismodes=1:17;
    data.edtmodes=[1:8,10];
    data.selmodes=1:6;
else
    % define 'compact'-mode
    fprintf('starting in "compact" mode...\n');
    set(handles.navppimgmode,'String',...
        {'normal','------','slices','labels','------','add#1','ind#1'});
    data.dismodes=[1,5,6,8,12,13,16];
    
    set(handles.edtppmode,'String',...
        {'none','major','------','proof'});
    data.edtmodes=[1,6,7,8];
    
    set(handles.navppselect,'String',{'none','normal','cluster','filter'});
    data.selmodes=[1,2,3,6];
end
% list of 'confirmed' conclusions
data.clist=[1,3,5,7];
% list of 'not sure' conclusions
data.nslist=[2,3,6,7];
% list of 'other' conslusions
data.otlist=[4,5,6,7];


% proof-notes
data.edfirst=[];        % first-slice-id of stack where corrections started
data.edlogid=[0,0];     % this holds [obj#,slice#] of the first correction
data.edconcl=0;         % this holds operator conclusion
data.edstats=[];        % this may hold additional stats, nothing now

% backup
data.internal=[];           % data for GUI-interface internal state
data.backup=25;             % how many proof-notes before "next" backup
data.bkcount=data.backup-1; % current counter for corrections
                            % set such that save is made before first cor.
if(flMonitor) data.bkcount=0; end                            
data.bksaves=11;            % number of backup saves so far, 
                            % make save full backup first time

% miscellaneous variables
data.fig1=[];           % handler of console window
data.fig2=[];           % handler for proof console window
data.cslc=1;            % slice # being edited
data.saved=0;           % whether backup had been saved
data.play=0;            % whether in playback mode
data.flip=0;            % whether in flip mode
data.curfirst=[];       % first section of the substack edited
data.undo=[];           % undo buffer
data.cmap=[];           % label->color map used
data.vXupdt=1;          % 'list of major processes' needs to be updated?
data.iproof=0;          % currently active proofing volume
data.ctime=zeros(1,7);  % time-stats array
data.axmainsize=[];     % size of graphic canvas for resizing
data.uislc=[];          % current slice for by-slice-filter modes
data.shadow=0;          % shadow mode flag
data.altcol=0;          % alternative color flag
data.modefirst=[];      % stack where proof list is initialized
data.K=1;               % size-reduction factor for display
data.overhead=3e7;        % overhead for memory control

% initialize counter ## of times move had been pressed, for ccdouble
%if(data.pref.ccdouble) data.cmove=0; else data.cmove=2; end
data.cmove=0;

data.seppos=get(handles.fmain,'Position'); % ref position for controls
end

% modify data.mmax for possible index variables in debug
if(isfield(debug,'ind1') && ~isempty(debug.ind1))
    for i=1:length(debug.ind1)
        if(~isempty(debug.ind1{i}))
            data.mmax=max(data.mmax,max(debug.ind1{i}(:)));
        end
    end
end

if(isfield(debug,'ind2') && ~isempty(debug.ind2))
    for i=1:length(debug.ind2)
        if(~isempty(debug.ind2{i}))
            data.mmax=max(data.mmax,max(debug.ind2{i}(:)));
        end
    end
end
data.mmax=double(data.mmax);

% initialize tmp variables
shape=[size(al{1}),length(al)];
img=zeros([shape(1:2),3],'uint8');      % primary image storage
tmp3=zeros([shape(1:2),3],'uint8');     % tmp image storage
idx=false(shape(1:2));                  % primary selection storage
tmp=zeros(shape(1:2),'uint32');         % multipurpose


% SET VARIABLE PROOF
% choose proof version
if(isempty(proof) || ~isfield(proof,'pmap'))
    proof.v=7;
elseif(~isfield(proof,'v'))
    fprintf('proof file version=3 or lower, assuming v=3...\n');    
    proof.v=3;
end

% load added corrections
if(isfield(proof,'first') && ~isempty(proof.first))
    s=['gui.buff',sprintf('%.3i',proof.first),'.mat'];
    if(firstuse && exist(s) && (isempty(flNoBuffLoad) || flNoBuffLoad))
        wcat=[];  data.tmp{20}=cell(size(cat));
        
        fprintf('found buffer file, loading\n');
        load(s);

        % apply draw-alterations stored in buff-file 'ind'
        for k=1:length(ind)
            if(~isempty(ind{k}))
                data.tmp{20}{k}=ind{k};
                if(~isempty(cat{k})) cat{k}(ind{k}{1})=ind{k}{2}; end
            end
        end
        
        if(~isempty(wcat)) debug.ind1=wcat; end
    end
end

% prepare proof.pmap -- cat0->cat assignments, import from debug if any
if(~isfield(proof,'pmap') || isempty(proof.pmap))
    if(~isempty(debug) && isfield(debug,'map'))
        proof.pmap=debug.map; 
    elseif(~isempty(debug) && isfield(debug,'mapping'))
        % THIS IS FOR BACKWARD COMPATIBILITY
        proof.pmap=debug.mapping;
    else        
        proof.pmap=uint32(linspace(0,data.mmax,data.mmax+1));
    end
elseif(length(proof.pmap)<data.mmax+1)
    proof.pmap(end+1:data.mmax+1)=length(proof.pmap):data.mmax;
end

% prepare proof.tmap -- objects id's
if(~isfield(proof,'tmap') || isempty(proof.tmap))
    proof.tmap=ones(size(proof.pmap),'uint8');
elseif(length(proof.tmap)<length(proof.pmap))
    a=ones(size(proof.pmap),'uint8');
    a(1:length(proof.tmap))=proof.tmap;
    proof.tmap=a;
end

% prepare proof.ttime -- proofing time & cor. statistics
% convention:
% [time between 'cycle', time for editing, classification id,...
%  navigation clicks, link clicks, quick edit clicks, draw/del clicks]
if(~isfield(proof,'ttime') || isempty(proof.ttime))
    proof.ttime=zeros(length(proof.pmap),7);
elseif(size(proof.ttime,1)<length(proof.pmap))
    a=zeros(length(proof.pmap),7);
    a(1:length(proof.ttime),:)=proof.ttime;
    proof.ttime=a;
end
    
if(proof.v<6)
    % THIS IS FOR BACKWARD COMPATIBILITY
    %   convention <6 uses n/a,dndr,spn,axn,glia,other
    %   convention >=6 uses extended types
    
    % check if conversion had previously occured
    if(isfield(proof,'note')) str=proof.note; else str={}; end
    if(~iscell(str)) str={str}; end

    % see if despite v<6, types had been converted before
    X=strfind(str,'itypes converted'); flg=0;
    for i=1:length(X) flg=flg | (~isempty(X{i})); end
    flg=flg & (max(proof.tmap)<=6);

    if(~flg)        
        % if here - means need to convert types:
        stats=regionprops(proof.tmap,'Area','PixelIdxList');
        for i=(length(stats)+1):6 stats(i).PixelIdxList=[]; end;
        
        proof.tmap(:)=1;      
        proof.tmap(stats(2).PixelIdxList)=2;        
        proof.tmap(stats(3).PixelIdxList)=8;
        proof.tmap(stats(4).PixelIdxList)=3;
        proof.tmap(stats(5).PixelIdxList)=4;
        proof.tmap(stats(6).PixelIdxList)=11;

        % add note about conversion
        str{length(str)+1}='itypes converted vX->v6';
        proof.note=str;
    end
end

if(~isfield(proof,'notes') || isempty(proof.notes))
    % prepare proof.notes
    proof.notes=cell(maxnote,5);
    proof.notes=cell2struct(proof.notes,{'tag','slice','concl','stats','note'},2);
    data.ednotes=0;
else
    tmp=proof.notes;    

    % THIS IS FOR BACKWARD COMPATIBILITY
    %   convention <7 uses different format for proof-notes    
    if(proof.v<7)         
        if(~isfield(debug,'logs'))
            fprintf('Warning: GUI v7 is incompatible with proofs produced\n');
            fprintf(' by earlier GUI versions; v<7 proof has been found!\n');
            fprintf(' Cannot convert vX->v7 because cannot find debug.logs!\n');
            fprintf(' Proof notes will be overwritten!\n');
            fprintf(' Supply debug.logs to convert proof vX->v7.\n');
            data.ednotes=0;
        else
            a=[];
            for i=1:length(tmp) if(~isempty(tmp(i).logid)) a=[a,i]; end; end
            proof.notes=cell(length(a)+maxnote,5);
            proof.notes=cell2struct(proof.notes,...
                {'tag','slice','concl','stats','note'},2);

            for i=1:length(a)
                % v<7 reference proof-notes by debug.logs events
                log=debug.logs(tmp(a(i)).logid);
                
                % obtain what is the type of reference event
                etype=log.event;
                
                % identify object # referenced in the event
                if(ismember(etype,[6,7])) nclst=log.tag;
                else nclst=log.outgoing; end
                nclst=proof.pmap(nclst+1);
                note.tag=nclst;
                
                % identify slice # referenced in the event
                note.slice=log.slice;
                % v<=3 uses different "conclusion" convention
                if(proof.v<=3)
                    switch(tmp(a(i)).concl)
                        case 2 % 2 means "confirmed"
                            note.concl=1;
                            note.stats=[];
                            note.note='';
                        case 3 % 3 means "not sure"
                            note.concl=2;
                            note.stats=[];
                            note.note=tmp(a(i)).note;
                        otherwise % otherwise assign "other"
                            note.concl=4;
                            note.stats=[];
                            note.note=tmp(a(i)).note;
                    end
                else
                    note.concl=tmp(a(i)).concl;
                    note.stats=tmp(a(i)).stats;
                    note.note=tmp(a(i)).note;
                end

                proof.notes(i)=note;
            end
            data.ednotes=length(a);            
        end
    else    
        % locate last written note
        a=[];
        for i=1:length(tmp) if(~isempty(tmp(i).tag)) a=[a,i]; end; end
        proof.notes=cell(length(a)+maxnote,5);
        proof.notes=cell2struct(proof.notes,...
                            {'tag','slice','concl','stats','note'},2);
        proof.notes(1:length(a))=tmp(a);
        data.ednotes=length(a);
    end
end
    
% assign first slice's global number
if(~isfield(proof,'first')) proof.first=1; end
if(~isfield(proof,'nstack')) proof.nstack=length(cat); end

% update massive length array
data.mmax=max(data.mmax,length(proof.pmap)-1);
data.mmax=max(data.mmax,double(max(proof.pmap)));

% assign proof-file version
data.version=proof.v;



% SET UP DEBUGING/PROOFING PROCESS
% parameters for computer-guided lost+found editing
% data.lfThrs defines which lost+found fragments are to be shown:
% [min anchors %% in the projection into neighbour section, 
%  min anchors %% of a member of the process in any one slice, 
%  min MinorAxis of a member of the process in any one slice, 
%  min length of process in ## of sections,
%  keep major processes which are at least this long in ## of sections,
%  min MajorAxis of the process: overrides anything else (not shown),
%  min MajorAxis that overrides everything else (shown)]
% major process is a process which reaches out to the volume boundary
data.lfThrs=[0.1,0.6,25,3,2,5,38];

% Processing debuging information
fdebug=0;           % flag whether debug had been updated
% assign .map field if in older version of debug (v3)
if(~isfield(debug,'map') & isfield(debug,'mapping'))
    fdebug=1;
    debug.map=debug.mapping;
end

if(isfield(debug,'logs'))
    fdebug=~isfield(debug,'idd') | ~isfield(debug,'idk') | ...
          ~isfield(debug,'idt') | ~isfield(debug,'ide');
      
    % converting to cells      
    if(isfield(debug,'slc') && ~iscell(debug.slc))
        fdebug=1;
        fprintf('GUI v7 uses image data in cell format, attempting to convert...\n');
        celltmp=cell(1,size(debug.slc,3));
        for k=1:length(celltmp) celltmp{k}=debug.slc(:,:,k); end
        debug.slc=celltmp; clear celltmp;
        fprintf('convert debug.slc successful\n');
    end
    if(isfield(debug,'cat0') && ~iscell(debug.cat0))
        fdebug=1;
        fprintf('GUI v7 uses image data in cell format, attempting to convert...\n');
        celltmp=cell(1,size(debug.cat0,3));
        for k=1:length(celltmp) celltmp{k}=debug.cat0(:,:,k); end
        debug.cat0=celltmp; clear celltmp;
        fprintf('convert debug.cat0 successful\n');
    end
    if(isfield(debug,'anchors') && ~iscell(debug.anchors))
        fdebug=1;
        fprintf('GUI v7 uses image data in cell format, attempting to convert...\n');
        celltmp=cell(1,size(debug.anchors,3));
        for k=1:length(celltmp) celltmp{k}=debug.anchors(:,:,k); end
        debug.anchors=celltmp; clear celltmp;
        fprintf('convert debug.anchors successful\n');
    end
    if(isfield(debug,'dbgimg') && ~iscell(debug.dbgimg))
        fdebug=1;
        fprintf('GUI v7 uses image data in cell format, attempting to convert...\n');
        celltmp=cell(1,size(debug.dbgimg,3));
        for k=1:length(celltmp) celltmp{k}=debug.dbgimg(:,:,k); end
        debug.dbgimg=celltmp; clear celltmp;
        fprintf('convert debug.dbgimg successful\n');
    end
    if(isfield(debug,'mskE') && ~iscell(debug.mskE))
        fdebug=1;
        fprintf('GUI v7 uses image data in cell format, attempting to convert...\n');
        celltmp=cell(1,size(debug.mskE,3));
        for k=1:length(celltmp) celltmp{k}=debug.mskE(:,:,k); end
        debug.mskE=celltmp; clear celltmp;
        fprintf('convert debug.dbgimg successful\n');
    end

    
    if(isfield(debug,'add1') && ~iscell(debug.add1))
        fprintf('GUI v7 uses image data in cell format...\n');
        fprintf('debug.add1 is not in cell format, clearing\n');
        debug.add1={};
    end
    if(isfield(debug,'add2') && ~iscell(debug.add2))
        fprintf('GUI v7 uses image data in cell format...\n');
        fprintf('debug.add2 is not in cell format, clearing\n');
        debug.add2={};
    end
    if(isfield(debug,'add3') && ~iscell(debug.add3))
        fprintf('GUI v7 uses image data in cell format...\n');
        fprintf('debug.add3 is not in cell format, clearing\n');
        debug.add3={};
    end
    if(isfield(debug,'ind1') && ~iscell(debug.ind1))
        fprintf('GUI v7 uses image data in cell format...\n');
        fprintf('debug.add1 is not in cell format, clearing\n');
        debug.ind1={};
    end
    if(isfield(debug,'ind2') && ~iscell(debug.ind2))
        fprintf('GUI v7 uses image data in cell format...\n');
        fprintf('debug.add1 is not in cell format, clearing\n');
        debug.ind2={};
    end            

    % check for missing fields:    
    % assign .anchors field if older version of debug (v3)
    if(~isfield(debug,'anchors') & force_fields)
        fdebug=1;
        fprintf('no debug.anchors found, building...\n');
        debug.anchors=cell(size(al));
        for k=1:data.smax debug.anchors{k}=im2bw(al{k},data.bgthr); end
    end
    
    % assign edges map
    if(~isfield(debug,'dbgimg') & force_fields)
        fdebug=1;
        fprintf('no debug.dbgimg found, building...\n');
        debug.dbgimg=cell(size(cat));
        for k=1:data.smax debug.dbgimg{k}=im2uint8(cat{k}==0); end
    end    

    % assign .idd field if older version of debug
    if(~isfield(debug,'idd'))
        fprintf('no debug.idd found, building...\n');
        debug.idd=debug.map(1+double([debug.logs.outgoing]));    
    end
    
    % assign .idk field if older version of debug
    if(~isfield(debug,'idk'))
        fprintf('no debug.idk found, building...\n');        
        debug.idk=[debug.logs.slice];
    end
    
    % assign .idt field if older version of debug
    if(~isfield(debug,'idt'))
        fprintf('no debug.idt found, building... \n');
        debug.idt=debug.map(1+double([debug.logs.tag]));
    end      
    
    % assign .ide field if older version of debug
    if(~isfield(debug,'ide'))
        fprintf('no debug.ide found, building... \n');
        debug.ide=[debug.logs.event];
    end
end
    
% Significant objects for guided lost+found search
% read additional info in stats for 3D objects in link.m
if(data.mmax==1) debug.unique=1; end
if(~isfield(debug,'unique'))
    fprintf('building list of unique 3D objects.');
    dn=ceil(data.smax/30);    
    
    fdebug=1;
    tmp=[];
    for i=1:length(cat) 
        if(rem(i,dn)==0) fprintf('.'); end
        tmp=union(tmp,unique(cat{i})); 
    end
    data.tmp{19}=setdiff(tmp,0);
    
    debug.unique=data.tmp{19};
    fprintf('\n');
else
    data.tmp{19}=debug.unique;
end

if(data.mmax==1) debug.major=1; end
if(~isfield(debug,'major'))
    fdebug=1;
    fprintf('no debug.major found, building.\n');   
    if(isfield(debug,'sstats') && ~isempty(debug.sstats))    
        tmp=debug.sstats;

        % category IIa: clusters that are too wide somewhere
        % i.e. MinorAxisLength is large somewhere
        ids=tmp(:,4)' > data.lfThrs(3);

        % category IIb: clusters that are too long in ## of sections
        ids=ids | (tmp(:,6)' > data.lfThrs(4));

        % category IIc: necessarily must be quite bright somewhere 
        % [%% anchors]
        ids=ids | (tmp(:,2)' > data.lfThrs(2));

        % category IIc: override: but not clusters that are too small
        ids=ids & (tmp(:,3)' > data.lfThrs(6));
        
        % category IId: override: if clusters are too large
        ids=ids | (tmp(:,3)' > data.lfThrs(7));


        % category III: HANDLING ESCAPING FRAGMENTS, MODIFY HERE
        % exclude major clusters
        %  with exception for those longer than ## slices
        if(isfield(debug,'esc'))
            tmp=debug.esc(tmp(debug.esc,6) < data.lfThrs(5));
        else
            tmp=[];
        end
        data.major=setdiff(find(ids),tmp);
        
        % only ids actually present in the stack somewhere
        data.major=intersect(data.major,data.tmp{19});
        data.major=data.major(:)';
        debug.major=data.major;
    else
        data.major=unique(proof.pmap(data.tmp{19}+1));
        data.major=data.major(:)';
        debug.major=data.major;
    end    
else
    data.major=debug.major(:)';
end

% define section to be shown for particular element of major array
if(~isfield(debug,'zmajor'))
    data.zmajor=zeros(size(data.major),'uint16');
else
    data.zmajor=debug.zmajor(:)';
end

% define ordering of major processes
if(data.mmax==1) debug.morder=1; end
if(~isfield(debug,'morder'))    
    if(prooforder)
        fdebug=1;
        fprintf('no debug.morder found, ordering proofing list.');
        dn=ceil(data.smax/30);
        morder=zeros(data.mmax,1);

        for i=1:data.smax
            if(rem(i,dn)==0) fprintf('.'); end
            tmp=cat{i};
            tmp(~ismember(proof.pmap(tmp+1),data.major))=0;

            tmp=sort(proof.pmap(tmp(:)+1));
            tmp=tmp(:);
            dtmp=diff([double(tmp);double(max(tmp))+1]);
            sums=diff(find([1;dtmp]));
            ind=tmp(dtmp>0);
            sums=sums(ind>0);
            ind=ind(ind>0);

            morder(ind)=morder(ind)+sums;
        end
        fprintf('\n');

        data.morder=morder(data.major);
        data.morder=data.morder(:)';
        debug.morder=data.morder;
    else
        % if prooforder==0, don't need to sort major's
        debug.morder=1:length(data.major);
    end
else
    data.morder=debug.morder(:)';
end

% create list of objects z-bounds
if(data.mmax==1) debug.objk=[1,1]; end
if(~isfield(debug,'objk'))
    fdebug=1;
    fprintf('no debug.objk found, building list of z-bounds.');
    kbounds=zeros(data.mmax,2);
    kbounds(:,1)=data.smax+1;
    dn=ceil(data.smax/30);    
    for i=1:data.smax
        if(rem(i,dn)==0) fprintf('.'); end
        ind=unique(cat{i});
        ind=setdiff(ind,0);
        kbounds(ind,1)=min(kbounds(ind,1),i);
        kbounds(ind,2)=max(kbounds(ind,2),i);
    end
    fprintf('\n');
    data.objk=kbounds;
    debug.objk=kbounds;
else
    data.objk=debug.objk;
end 

% volume information for significant objects
if(~isfield(debug,'vorder')) 
    data.vorder=data.morder; 
else
    data.vorder=debug.vorder; 
end

% modify ordering list to show first longer objects
global zRelOrder
if(isempty(zRelOrder)) zRelOrder=1000; end
v=max(0,data.objk(data.major,2)-data.objk(data.major,1)-2)';
data.morder=v*zRelOrder+data.morder;


% list of escaping objects
if(~isfield(debug,'esc')) fdebug=1; debug.esc=[]; end

% print warning that debug had been updated
if(fdebug)
    if(firstuse)
        fprintf('debug variable had been modified, consider update-saving\n');
        fprintf('debug variable in the primary data file\n');
        fprintf('use save filename -append debug\n');
    else
        fprintf('updating debug variable in the data file\n');
        save(eventdata,'debug','-append');
    end
end
    
% data.smap holds slice ## where the process was last updated
data.smap=uint16(zeros(size(proof.pmap)));

tic;
if(firstuse) fprintf('initialization is complete\n'); end
% -----------------------------------------------------------------------



% -----------------------------------------------------------------------
% set slider browser scale
set(handles.navsldbrowse,'Min',1,'Max',data.smax);
set(handles.navsldbrowse,'SliderStep',[1/(data.smax-1),...
                                        1/max(10,(data.smax-1)/10)]);
% reset slider browser only if at first use
if(firstuse) set(handles.navsldbrowse,'Value',1); end

% set info-display for the first section
set(handles.navtxtslice,'String',sprintf('Section # %i*',...
    proof.first+get(handles.navsldbrowse,'Value')-1));                 

% see if there is backup
% COMPATIBILITY SWITCH
s_old='gui.backup.mat';
if(exist('gui.backup.mat'))
    tmp=proof;
    load gui.backup.mat proof
    if(tmp.first==proof.first)
        s=sprintf('gui.backup%.3i.mat',proof.first);
        movefile('gui.backup.mat',s);
    else
        proof=tmp;
    end
end
% COMPATIBILITY SWITCH        

s=sprintf('gui.backup%.3i.mat',proof.first);
if(exist(s) && firstuse)
    % check that backup is for currently loaded substack
    tmp=proof;
    load(s,'proof');
    if(tmp.first==proof.first)
        x=questdlg('Found gui.backup.mat, load old session?');               
        % x=input('found gui.backup.mat, load old session [y/n]?','s');
        % if(strcmp(x,'y') | strcmp(x,'Y')) savebackup(1,handles); end
        if(strcmp(x,'Yes')) savebackup(1,handles); else proof=tmp; end
        if(strcmp(x,'Cancel')) error('User requested termination!'); end        
    else
        proof=tmp;
    end
end
                                        
% DISPLAY
utdgetimg(handles,[]);
utdgetidx(handles,[]);
utdshow(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% =======================================================================
%                       AUXILIARY UTILITIES
% =======================================================================

% #######################################################################
% SOMETHING MAY GO HERE




% =======================================================================
%                           DRAWING AND NAVIGATION

% #######################################################################
function utdgetimg(handles,box)
% this function is responsible for building the image of the section drawn;
%  the image either uint8(:,:,3) true color array for color or uint8(:,:)
%  for bw is stored in global variable img

% notes:
% may change storage to uint32(:,:) indexed array if figure out how to
%  draw uint32 indexed color images???

% global data storages
global al cat debug proof

% shared tmp variables
global tmp img idx tmp3

% shared internal data
global data
% initializing
k=round(get(handles.navsldbrowse,'Value')); % current section id
ifget('alcat',k);
ik=k-data.cpos+1;
p=get(handles.navsldmixer,'Value');         % current mixer value
zoom=data.zoom;                             % current display window
shade=data.pref.shade(1);                        % set shading level

imap=get(handles.navrbtmap,'Value');        % is there additional mapping?
ioverlap=get(handles.navckoverlap,'Value'); % are we showing overlap?
ioverlay=get(handles.navckoverlay,'Value'); % are we showing overlay?
iflip=data.flip;                            % are we in flip mode?

mshow=get(handles.navppimgmode,'Value');    
mshow=data.dismodes(mshow);                 % what is the display mode?
mselect=get(handles.navppselect,'Value');   
mselect=data.selmodes(mselect);             % what is the selection mode?
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);                   % is there editing-selection?

itrack=(mselect>1);                         % are we tracking?

ids=data.ids;                               % selection clusters

logid=data.logid;                           % active log record
if(logid>0)
    etype=debug.logs(logid).event;          % event id, if any
    kk=debug.logs(logid).slice;             % event section
else
    etype=-1;                               % otherwise non 
    kk=-1;
end

% set drawing mode for first overlay
if((ioverlap)&(k>1)&(~ioverlay))
    ovk=1;
elseif(ioverlap & ~ioverlay)
    ovk=2;
else
    ovk=3;
end

% -----------------------------------------------------------------------
% define operation box
if(nargin<2)            % if none specified - current zoom
    box={zoom(1):zoom(2),zoom(3):zoom(4)};
elseif(isempty(box))    % if [] specified - full zoom
    box={1:data.shape(1),1:data.shape(2)};
end
% set processed zoom region
data.v1zoom=[box{1}(1),box{1}(end),box{2}(1),box{2}(end)];
lshape=size(cat{k}(box{:}));

% are we in flip-mode
% EXPERIMENTAL SHADOW ADDON HERE
if(data.shadow & ~iflip)
    if(strcmp(get(handles.uimshadowmbr,'Checked'),'off'))
        if(length(size(img))>2 || size(img,3)>1) 
            img=zeros(data.shape,'uint8'); 
        end
        img(box{:})=immultiply(al{ik}(box{:}),...
            1/max(data.pref.shade));
        return
%         p=1; mshow=6; ioverlay=0; ovk=3;
    else    
        if(size(img,3)<3) img=zeros([data.shape,3],'uint8'); end
        tmp=immultiply(al{ik}(box{:}),1/max(data.pref.shade));
        for i=1:3 img(box{:},i)=tmp; end
        % ind2 will contain 1pxl separators
        if(isfield(debug,'ind2') && ~isempty(debug.ind2) && ...
                ~isempty(debug.ind2{ik}))
            tmp(debug.ind2{ik}(box{:})==0)=...
                imadd(tmp(debug.ind2{ik}(box{:})==0),double(im2uint8(p)));

        else
            tmp(cat{ik}(box{:})==0)=...
                imadd(tmp(cat{ik}(box{:})==0),double(im2uint8(p)));
        end
        img(box{:},2)=tmp;
        return;
    end    
end
    
if(iflip) p=1; mshow=6; ioverlay=0; ovk=3; shadow=1; itrack=0; end

% First layer image
tmp=[];
switch(ovk)
    case 1      % drawing bw overlap
        switch(mshow)
            case {1,7,8}  % regular/labels
                tmp=imlincomb(p,cat{ik}(box{:})>0,1-p,cat{ik-1}(box{:})>0);
            case 6      % min-overlay of slices
                tmp=min(al{ik}(box{:}),al{ik-1}(box{:}));
            case {2,9}  % regular/forward
                if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=imlincomb(p,debug.cat0{k}(box{:})>0,...
                                        1-p,debug.cat0{k-1}(box{:})>0);
                end
            case {3,10} % regular anchors
                if(isfield(debug,'anchors') && ~isempty(debug.anchors))
                    tmp=imlincomb(p,debug.anchors{k}(box{:})>0,...
                                        1-p,debug.anchors{k-1}(box{:})>0);
                end
            case {4,11} % regular/membranes
                if(isfield(debug,'dbgimg') && ~isempty(debug.dbgimg))
                    tmp=imlincomb(p,debug.dbgimg{k}(box{:})>0,...
                                    1-p,debug.dbgimg{k-1}(box{:})>0);
                end
            case 13     % add1
                if(isfield(debug,'add1') && ~isempty(debug.add1)...
                        && ~isempty(debug.add1{k}) && ~isempty(debug.add1{k-1}))
                    tmp=imlincomb(p,debug.add1{k}(box{:})>0,...
                                     1-p,debug.add1{k-1}(box{:})>0);
                end
            case 14     % add2
                if(isfield(debug,'add2') && ~isempty(debug.add2)...
                        && ~isempty(debug.add2{k}) && ~isempty(debug.add2{k-1}))
                    tmp=imlincomb(p,debug.add2{k}(box{:})>0,...
                                     1-p,debug.add2{k-1}(box{:})>0);
                end
            case 15     % add3
                if(isfield(debug,'add3') && ~isempty(debug.add3)...
                        && ~isempty(debug.add3{k}) && ~isempty(debug.add3{k-1}))
                    tmp=imlincomb(p,debug.add3{k}(box{:})>0,...
                                    1-p,debug.add3{k-1}(box{:})>0);
                end
            case 16     % add4
                ifget('wcat',k);
                if(isfield(debug,'ind1') && ~isempty(debug.ind1)...
                        && ~isempty(debug.ind1{k}) && ~isempty(debug.ind1{k-1}))
                    tmp=imlincomb(p,debug.ind1{k}(box{:})>0,...
                        1-p,debug.ind1{k-1}(box{:})>0);
                end
            case 17     % add5
                if(isfield(debug,'ind2') && ~isempty(debug.ind2)...
                        && ~isempty(debug.ind2{k}) && ~isempty(debug.ind2{k-1}))
                    tmp=imlincomb(p,debug.ind2{k}(box{:})>0,...
                        1-p,debug.ind2{k-1}(box{:})>0);
                end
        end
    case 2          % drawing overlap in an edge secton
        switch(mshow)
            case {1,7,8}
                tmp=cat{ik}(box{:})>0;
            case 6      % slices
                tmp=al{ik}(box{:});
            case {2,9}
                if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=debug.cat0{k}(box{:})>0;
                end
            case {3,10}
                if(isfield(debug,'anchors') && ~isempty(debug.anchors))
                    tmp=debug.anchors{k}(box{:})>0;
                end
            case {4,11}
                if(isfield(debug,'dbgimg') && ~isempty(debug.dbgimg))
                    tmp=debug.dbgimg{k}(box{:})>0;
                end
            case 13
                if(isfield(debug,'add1') && ~isempty(debug.add1)...
                        && ~isempty(debug.add1{k}))
                    tmp=debug.add1{k}(box{:})>0;
                end
            case 14
                if(isfield(debug,'add2') && ~isempty(debug.add2)...
                        && ~isempty(debug.add2{k}))
                    tmp=debug.add2{k}(box{:})>0;
                end
            case 15
                if(isfield(debug,'add3') && ~isempty(debug.add3)...
                        && ~isempty(debug.add3{k}))
                    tmp=debug.add3{k}(box{:})>0;
                end
            case 16
                ifget('wcat',k);                
                if(isfield(debug,'ind1') && ~isempty(debug.ind1)  && ~isempty(debug.ind1{k}))
                    tmp=debug.ind1{k}(box{:})>0;
                end
            case 17
                if(isfield(debug,'ind2') && ~isempty(debug.ind2)  && ~isempty(debug.ind2{k}))
                    tmp=debug.ind2{k}(box{:})>0;
                end
        end
    case 3      % drawing datasets themselves
        switch(mshow)
            case {1,8}
                tmp=cat{ik}(box{:});
            case {2,9}
                if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=debug.cat0{k}(box{:});
                end
            case {3,10}
                if(isfield(debug,'anchors') && ~isempty(debug.anchors))
                    tmp=debug.anchors{k}(box{:});
                end
            case {4,11}
                if(isfield(debug,'dbgimg') && ~isempty(debug.dbgimg))
                    tmp=debug.dbgimg{k}(box{:});
                end
            case 6      % slices
                tmp=al{ik}(box{:});
            case 7
                tmp=cat{ik}(box{:})>0;
            case 13
                if(isfield(debug,'add1') && ~isempty(debug.add1)...
                        && ~isempty(debug.add1{k}))
                    tmp=debug.add1{k}(box{:});
                end
            case 14
                if(isfield(debug,'add2') && ~isempty(debug.add2)...
                        && ~isempty(debug.add2{k}))
                    tmp=debug.add2{k}(box{:});
                end
            case 15
                if(isfield(debug,'add3') && ~isempty(debug.add3)...
                        && ~isempty(debug.add3{k}))
                    tmp=debug.add3{k}(box{:});
                end
            case 16
                ifget('wcat',k);                
                if(isfield(debug,'ind1') && ~isempty(debug.ind1) && ~isempty(debug.ind1{k}))
                    tmp=debug.ind1{k}(box{:});
                end
            case 17
                if(isfield(debug,'ind2') && ~isempty(debug.ind2) && ~isempty(debug.ind2{k}))
                    tmp=debug.ind2{k}(box{:});
                end
        end
end
% if caught nothing, assign zeros
if(isempty(tmp)) tmp=zeros(lshape,'uint8'); end


% see if we are drawing a color image
imgsize=ismember(mshow,[1,2,8,9,16,17]) & (ovk==3);
% if color image, store colors as our first layer
altcol=strcmp(get(handles.uimmult,'Checked'),'on');
if(imgsize)
    if(isempty(data.cmap) || data.mmax>size(data.cmap,1)-1 || ...
            max(proof.pmap)>size(data.cmap,1)-1)
        zerocolor=[0 0 0];
        if(isstr(data.pref.clmp)) 
            cmap=feval(data.pref.clmp,data.mmax); 
        else cmap=data.pref.clmp; end
        % shuffle
        S = rand('state');
        rand('state', 0);
        index = randperm(data.mmax);
        cmap = cmap(index,:,:);
        rand('state', S);
        % form cmap
        data.cmap = im2uint8([zerocolor;cmap]);
    end
%     % apply color-remaping if requested    
%     if(imap) tmp=proof.pmap(imadd(tmp,1)); end
    clmp=data.cmap;
    if(imap) clmp=clmp(proof.pmap(1:length(proof.pmap))+1,:); end
    if(altcol)
        clmp=im2double(clmp);  clmp(1,:)=[1 1 1];

        % adjust colors      
%         x=sum(double(clmp>0),2);
%         x(x==1)=0.5; x(x==2)=0.75; x(x==3)=1;
%         x=x./max(clmp,[],2);
%         x=x./sum(clmp,2);
        x=1./max(clmp,[],2);

        clmp=clmp.*repmat(x,[1,3]);
%        clmp=clmp./repmat(sum(clmp.^1,2),[1,3]);        
        clmp(1,:)=[1 1 1]*p; 
        clmp=im2uint8(clmp);
    end
    
    tmp3=zeros([size(tmp),3],'uint8');
    
    %=========================================================
    % EXPERIMENTAL IMAGE REDUCTION MODE
    sm=get(handles.axmain,'Position');
    sm=sm([4,3])-sm([2,1]);
    K=max(1,min(length(box{1})/sm(1),length(box{2})/sm(2)));
    
    tmp1=imresize(tmp,1/K,'nearest');
    tmp2=zeros([size(tmp1),3],'uint8');
    tmp2(:)=[clmp(1+tmp1(:),1);clmp(1+tmp1(:),2);clmp(1+tmp1(:),3)];

    tmp2=imresize(tmp2,K,'nearest');
    sn=min(size(tmp2),size(tmp3));
    tmp3(1:sn(1),1:sn(2),:)=tmp2(1:sn(1),1:sn(2),:);
    % EXPERIMENTAL IMAGE REDUCTION MODE
    %==========================================================
    
%     tmp3(:)=[clmp(1+tmp(:),1);clmp(1+tmp(:),2);clmp(1+tmp(:),3)];
else
    % if not, just convert to grayscale
    tmp3=im2uint8(tmp);
end

% see if we are tracking, set shading levels
if(itrack) ishade=shade; else ishade=1; end
    
% Second layer image
if(~ioverlay & (ovk==3) & ismember(mshow,[1,2]))
    % if labels and grayscale
    if(size(img,3)<3) img=zeros([data.shape,3],'uint8'); end 

    if(altcol)
        for i=1:3
            tmp3(:,:,i)=im2uint8(im2double(tmp3(:,:,i)).*...
                                           im2double(al{ik}(box{:}))); 
        end
        x=1/ishade/im2double(max(tmp3(:)));
        tmp3(:)=immultiply(tmp3(:),x);
        for i=1:3 img(box{:},i)=tmp3(:,:,i); end
    else
        %=========================================================
        % EXPERIMENTAL IMAGE REDUCTION MODE
        tmp1=imresize(al{ik}(box{:}),1/K,'nearest');
        tmp1=repmat(tmp1,[1 1 3]);
        
        tmp2=imresize(tmp3,1/K,'nearest');
        tmp1=imlincomb(p/ishade,tmp1,(1-p)/ishade,tmp2);
        tmp1=imresize(tmp1,K,'nearest');
        tmp2=img(box{:},:);
        sn=min(size(tmp2),size(tmp1));        
        tmp2(1:sn(1),1:sn(2),:)=tmp1(1:sn(1),1:sn(2),:);
        img(box{:},:)=tmp2;        
        % EXPERIMENTAL IMAGE REDUCTION MODE
        %==========================================================
        
        
%         for i=1:3
%             img(box{:},i)=imlincomb(p/ishade,al{ik}(box{:}),...
%                                             (1-p)/ishade,tmp3(:,:,i));
%         end
    end
elseif(~ioverlay & (ovk==3) & ismember(mshow,[8,9,16,17]))
    % if any other labels color dataset
    if(size(img,3)<3) img=zeros([data.shape,3],'uint8'); end    
    
    for i=1:3  img(box{:},i)=immultiply(tmp3(:,:,i),1/ishade);   end
elseif(~ioverlay)
    % in any other case & not overlay show grayscale image
    if(size(img,3)>2) img=zeros(data.shape,'uint8'); end
    
    img(box{:})=immultiply(tmp3,1/ishade);
elseif(ioverlay & ismember(mshow,[1,2,3,4]))
    % special overlay: chklovskii ThreeColor overlay
    if(size(img,3)<3) img=zeros([data.shape,3],'uint8'); end    
    
    for i=1:3  
        if((k+i-2>0) & (k+i-3<data.smax))
            img(box{:},i)=imcomplement(al{ik-2+i}(box{:})); 
        else
            img(box{:},i)=im2uint8(zeros(lshape));  
        end
    end
elseif(ioverlay & ~isempty(data.tmp{1}))
    % if we do do overlay, define what we want to draw in second layer
    tmp=[];
    switch(data.tmp{1})
        case 8
            tmp=cat{ik}(box{:});
        case 9
            if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                tmp=debug.cat0{k}(box{:});
            end
        case 10
            if(isfield(debug,'anchors') && ~isempty(debug.anchors))
                tmp=debug.anchors{k}(box{:});
            end
        case 11
            if(isfield(debug,'dbgimg') && ~isempty(debug.dbgimg))
                tmp=debug.dbgimg{k}(box{:});
            end
        case 6
            tmp=al{ik}(box{:});
        case 7
            tmp=cat{ik}(box{:})>0;
        case 13
            if(isfield(debug,'add1') && ~isempty(debug.add1)...
                    && ~isempty(debug.add1{k}))
                tmp=debug.add1{k}(box{:});
            end
        case 14
            if(isfield(debug,'add2') && ~isempty(debug.add2)...
                    && ~isempty(debug.add2{k}))
                tmp=debug.add2{k}(box{:});
            end
        case 15
            if(isfield(debug,'add3') && ~isempty(debug.add3)...
                    && ~isempty(debug.add3{k}))
                tmp=debug.add3{k}(box{:});
            end
        case 16
            ifget('wcat',k);                            
            if(isfield(debug,'ind1') && ~isempty(debug.ind1)  && ~isempty(debug.ind1{k}))
                tmp=debug.ind1{k}(box{:});
            end
        case 17
            if(isfield(debug,'ind2') && ~isempty(debug.ind2)  && ~isempty(debug.ind2{k}))
                tmp=debug.ind2{k}(box{:});
            end
    end
    % if nothing, make zeros
    if(isempty(tmp)) tmp=zeros(lshape,'uint8'); end
    
    switch (data.tmp{1})
        case {8,9,16,17}    % drawing indexed colormap
            if(size(img,3)<3) 
                img=zeros([data.shape,3],'uint8'); 
            end
            
            % apply mapping if needed
            if(imap) tmp=proof.pmap(imadd(tmp,1)); end
            
            if(isempty(data.cmap))
                zerocolor=[0 0 0];
                if(isstr(data.pref.clmp)) 
                    cmap=feval(data.pref.clmp,data.mmax); 
                else
                    cmap=data.pref.clmp; 
                end
                % shuffle
                S = rand('state');
                rand('state', 0);
                index = randperm(data.mmax);
                cmap = cmap(index,:,:);
                rand('state', S);
                % form cmap
                data.cmap = im2uint8([zerocolor;cmap]);
            end
            img(box{:},:)=reshape([data.cmap(1+tmp(:),1);data.cmap(1+tmp(:),2);...
                    data.cmap(1+tmp(:),3)],[size(tmp,1),size(tmp,2),3]);
            
            if(imgsize)
                for i=1:3
                    img(box{:},i)=imlincomb(p/ishade,img(box{:},i),...
                                                (1-p)/ishade,tmp3(:,:,i));
                end
            else
                for i=1:3
                    img(box{:},i)=imlincomb(p/ishade,img(box{:},i),...
                                                        (1-p)/ishade,tmp3);
                end
            end                
        otherwise       % otherwise overlay with gray
            tmp=im2uint8(tmp);

            if(imgsize)
                if(size(img,3)<3)
                    img=zeros([data.shape,3],'uint8');
                end
                for i=1:3
                    img(box{:},i)=imlincomb(p/ishade,tmp,(1-p)/ishade,tmp3(:,:,i));
                end
            else
                if(size(img,3)>2) img=zeros(data.shape,'uint8'); end
                img(box{:})=imlincomb(p/ishade,tmp,(1-p)/ishade,tmp3);
            end
    end
end
        

% ######################################################################
function utdgetidx(handles,box);
% function is responsible for specifying shades selection, the shade
%   selection is stored in global variable idx as logical(:,:) array

% global data variables
global al cat debug proof

% shared tmp variables
global tmp img idx

% shared internal data
global data

% initializing
k=round(get(handles.navsldbrowse,'Value')); % current section id
ifget('cat',k);
ik=k-data.cpos+1;
p=get(handles.navsldmixer,'Value');         % current mixer value
zoom=data.zoom;                             % current display window
shade=data.pref.shade(1);                        % set shading level

imap=get(handles.navrbtmap,'Value');        % is there additional mapping?
ioverlap=get(handles.navckoverlap,'Value'); % are we showing overlap?
iflip=data.flip;      % are we in flip mode?

mshow=get(handles.navppimgmode,'Value');    
mshow=data.dismodes(mshow);                 % what is the display mode?
mselect=get(handles.navppselect,'Value');   
mselect=data.selmodes(mselect);             % what is the selection mode?
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);                   % is there editing-selection?

itrack=(mselect>1);                         % are we tracking?

ids=data.ids;                               % selection clusters

logid=data.logid;                           % current log record
if(logid>0)
    etype=debug.logs(logid).event;          % event id, if any
    kk=debug.logs(logid).slice;             % event section
    kk=kk-data.cpos+1;
else
    etype=-1;                               % otherwise zero 
    kk=-1;
end

% -----------------------------------------------------------------------
% define operation box
if((nargin<2))          % if none specified - current zoom
    box={zoom(1):zoom(2),zoom(3):zoom(4)};
elseif(isempty(box))    % if [] specified - full zoom
    box={1:data.shape(1),1:data.shape(2)};
end
% processed zoom range for selection
data.v2zoom=[box{1}(1),box{1}(end),box{2}(1),box{2}(end)];

% define flip-mode
if(iflip) p=1; mshow=1; end                


%=========================================================
% EXPERIMENTAL IMAGE REDUCTION MODE
sm=get(handles.axmain,'Position');
sm=sm([4,3])-sm([2,1]);
K=max(1,min(length(box{1})/sm(1),length(box{2})/sm(2)));
% EXPERIMENTAL IMAGE REDUCTION MODE
%==========================================================



% choose what to highlight
ff=0;
tmp=[];
switch(mselect)
    case 2      % normal selection mode, context sensitive
        switch(medt)
            case {1,6}      % select clusters
                %=========================================================
                % EXPERIMENTAL IMAGE REDUCTION MODE
                tmp1=cat{ik}(box{:});
                tmp=false(size(tmp1));
                
                tmp1=imresize(tmp1,1/K,'nearest');
                if(imap)
                    tmp1=ismember(proof.pmap(tmp1+1),proof.pmap(ids+1));
                else
                    tmp1=ismember(tmp1,ids);
                end                
                tmp1=imresize(tmp1,K,'nearest');
                
                sn=min(size(tmp),size(tmp1));
                tmp(1:sn(1),1:sn(2),:)=tmp1(1:sn(1),1:sn(2),:);
                % EXPERIMENTAL IMAGE REDUCTION MODE
                %==========================================================
                
                
%                 if(imap)    % use identification through proof.map if imap
%                     tmp=ismember(proof.pmap(imadd(cat{ik}(box{:}),1)),...
%                         proof.pmap(1+ids));
%                 else        % use original labels otherwise
%                     tmp=ismember(cat{ik}(box{:}),ids);
%                 end
            case 2  % select splits
                % select events
                idd=find(debug.idk==k);
                idd=idd(debug.ide(idd)==6);
                % extract tags
                idd=[debug.logs(idd).tag];
                % references here are immediately to cat0
                if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=ismember(debug.cat0{k}(box{:}),idd);
                else
                    tmp=ismember(cat{ik}(box{:}),idd);
                end
            case 3  % select mergers, see case 2 for explanations
                idd=find(debug.idk==k);
                idd=idd(debug.ide(idd)==3);
                % references here are to .slc or to debug.map'ed cat0
                if(isfield(debug,'slc') && ~isempty(debug.slc))
                    tmp=ismember(debug.slc{k}(box{:}),[debug.logs(idd).tag]);
                elseif(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=ismember(debug.map(debug.cat0{k}(box{:})+1),debug.idd(idd));
                else
                    tmp=ismember(debug.map(cat{ik}(box{:})+1),debug.idd(idd));
                end
            case 4  % select losts
                idd=find(debug.idk==k);
                idd=idd(debug.ide(idd)==7);
                idd=[debug.logs(idd).tag];
                if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=ismember(debug.cat0{k}(box{:}),idd);
                else
                    tmp=ismember(cat{ik}(box{:}),idd);
                end
            case 5  % select founds
                idd=find(debug.idk==k);
                idd=idd(debug.ide(idd)==9);
                idd=[debug.logs(idd).ccidx];
                if(isfield(debug,'slc') && ~isempty(debug.slc))
                    tmp=ismember(debug.slc{k}(box{:}),idd);
                elseif(isfield(debug,'cat0') && ~isempty(debug.cat0))
                    tmp=ismember(debug.cat0{k}(box{:}),debug.idd(idd));
                else
                    tmp=ismember(cat{ik}(box{:}),debug.idd(idd));
                end
%             case 6   % select major
%                 if(imap)    % use identification through proof.map if imap
%                     tmp=ismember(proof.pmap(cat{k}(box{:})+1),...
%                                                     proof.pmap(1+ids));
%                 else        % use original labels otherwise
%                     tmp=ismember(cat{ik}(box{:}),ids);
%                 end                
            case 8  % proof mode
                % extract major fragments
                if(data.vXupdt & get(handles.edtrbtredraw,'Value')) 
                    proof.major=proof.pmap(dbggetmajor+1);
                    proof.major=unique(proof.major);
                    data.vXupdt=0;
                end
                idd=proof.major;
                % select by reference to pmap-transformed cat
                tmp=ismember(proof.pmap(cat{ik}(box{:})+1),idd);
        end
    case 3 % 2 means select cluster no matter what, see above for imap
        %=========================================================
        % EXPERIMENTAL IMAGE REDUCTION MODE
        tmp1=cat{ik}(box{:});
        tmp=false(size(tmp1));

        tmp1=imresize(tmp1,1/K,'nearest');
        if(imap)
            tmp1=ismember(proof.pmap(tmp1+1),proof.pmap(ids+1));
        else
            tmp1=ismember(tmp1,ids);
        end
        tmp1=imresize(tmp1,K,'nearest');

        sn=min(size(tmp),size(tmp1));
        tmp(1:sn(1),1:sn(2),:)=tmp1(1:sn(1),1:sn(2),:);
        % EXPERIMENTAL IMAGE REDUCTION MODE
        %==========================================================
        
%         if(imap)
%             tmp=ismember(proof.pmap(cat{ik}(box{:})+1),proof.pmap(1+ids));
%         else
%             tmp=ismember(cat{ik}(box{:}),ids);
%         end
    case 4      % select participating members of active event, context sensitive
        switch(etype)
            case 6  % event is split
                if(kk<k)
                    % if we are above the event's slice, select split
                    % fragments by reference to slc
                    if(isfield(debug,'slc') && ~isempty(debug.slc))
                        tmp=ismember(debug.slc{kk+1}(box{:}),debug.logs(logid).ccidx);
                    elseif(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.map(debug.cat0{k}(box{:})+1),debug.idd(logid));
                    else
                        tmp=ismember(debug.map(cat{ik}(box{:})+1),debug.idd(logid));
                    end
                elseif(kk==k)
                    % if in the event's slice, highlight the main fragment 
                    % itself by reference to cat0
                    if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.cat0{kk}(box{:}),debug.logs(logid).tag);
                    else
                        ikk=min(data.csize,max(1,kk-data.cpos+1));
                        tmp=ismember(cat{ikk}(box{:}),debug.logs(logid).tag);
                    end
                else
                    % if below event's slice, highlight the main fragment
                    % there
                    if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.map(debug.cat0{k}(box{:})+1),debug.idt(logid));
                    else
                        tmp=ismember(debug.map(cat{ik}(box{:})+1),debug.idt(logid));
                    end
                end
            case 7 % event is lost
                if(kk<=k)
                    % if we are above or at the event's slice, select
                    % projection of main fragment
                    if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.cat0{kk}(box{:}),debug.logs(logid).tag);
                    else
                        ikk=min(data.csize,max(1,kk-data.cpos+1));
                        tmp=ismember(cat{ikk}(box{:}),debug.logs(logid).tag);
                    end
                else
                    % if below event's slice, highlight the main fragment
                    % there
                    if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.map(debug.cat0{k}(box{:})+1),debug.idt(logid));
                    else
                        tmp=ismember(debug.map(cat{ik}(box{:})+1),debug.idt(logid));
                    end
                end
            case 9  % event is found
                if(kk>=k)
                    % if we are below or at the event's slice, select
                    % projection of main fragment
                    if(isfield(debug,'slc') && ~isempty(debug.slc))
                        tmp=ismember(debug.slc{kk}(box{:}),debug.logs(logid).tag);
                    elseif(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.cat0{kk}(box{:}),debug.logs(logid).outgoing);
                    else
                        ikk=min(data.csize,max(1,kk-data.cpos+1));
                        tmp=ismember(cat{ikk}(box{:}),debug.logs(logid).outgoing);
                    end
                else
                    % if we are above the event's slice, select the main
                    % fragment there
                    if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.map(debug.cat0{k}(box{:})+1),debug.idd(logid));
                    else
                        tmp=ismember(debug.map(cat{ik}(box{:})+1),debug.idd(logid));
                    end
                end
            otherwise
                if(kk<=k && etype>0)
                    % if we are above the event's slice, select main
                    % fragment itself
                    if(isfield(debug,'slc') && ~isempty(debug.slc))
                        tmp=ismember(debug.slc{kk}(box{:}),debug.logs(logid).tag);
                    elseif(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.map(debug.cat0{k}(box{:})+1),debug.idd(logid));
                    else
                        tmp=ismember(debug.map(cat{ik}(box{:})+1),debug.idd(logid));
                    end
                elseif(etype>0)
                    % if we are below event's slice, highlight active 
                    % members of the event in the current slice by
                    % reference in ccidx->cat0
                    nflg=size(debug.logs(logid).flags,1);
                    flags=debug.logs(logid).flags(nflg,:);
                    ccids=debug.logs(logid).ccidx(flags);
                    if(isfield(debug,'cat0') && ~isempty(debug.cat0))
                        tmp=ismember(debug.cat0{k}(box{:}),ccids);
                    else
                        tmp=ismember(cat{ik}(box{:}),ccids);
                    end
                end
        end
    case 6      % 'filter' mode: filter for proof-read
        % need to have nonempty selection in data.tmp{2} and be in the
        % slice of the edited records
        if(~isfield(data,'ind') || isempty(data.tmp{7}))
            tmp=false(size(cat{1}(box{:})));
        else
            if(imap)
                tmp=ismember(proof.pmap(cat{ik}(box{:})+1),...
                                proof.pmap(data.tmp{7}(data.ind,1)+1));                
            else
                tmp=ismember(cat{ik}(box{:}),data.tmp{7}(data.ind,1));                
            end
        end
        
%         if(~isempty(data.tmp{2}) && ~isempty(data.tmp{3}) && (k==data.tmp{3}))
%             tmp=ismember(proof.pmap(cat{ik}(box{:})+1),data.tmp{2});
%         else
%             if(imap)  tmp=ismember(proof.pmap(imadd(cat{k}(box{:}),1)),proof.pmap(1+ids));
%                 else tmp=ismember(cat{ik}(box{:}),ids);  end
%         end
end
if(isempty(tmp)) tmp=false(size(cat{1}(box{:}))); end
% if selection was "same", need not update selection area
if(mselect~=5) idx(box{:})=tmp; end

% ######################################################################
function utdshow(handles)
% this function is responsible for actual drawing. It uses global variables
%   img and idx where image and shading are respectively stored

% global tmp variables
global img idx tmp tmp3

% global data variables
global data

% intialize stat
zoom=data.zoom;
v1zoom=data.v1zoom;
v2zoom=data.v2zoom;
% fov box
box={zoom(1):zoom(2),zoom(3):zoom(4)};

% if out of processed windows, need to redraw!!!
iv1=isempty(v1zoom);
if(~iv1)
    iv1=((zoom(1)<v1zoom(1)) | (zoom(2)>v1zoom(2)) |...
        (zoom(3)<v1zoom(3)) | (zoom(4)>v1zoom(4)));
end

iv2=isempty(v2zoom);
if(~iv2)
    iv2=((zoom(1)<v2zoom(1)) | (zoom(2)>v2zoom(2)) |...
        (zoom(3)<v2zoom(3)) | (zoom(4)>v2zoom(4)));
end

% initialization
mshow=get(handles.navppimgmode,'Value');    
mshow=data.dismodes(mshow);                 % what is the display mode?
iedge=(mshow==4);                           % are we showing edges?
iflip=data.flip;                            % are we in flip mode?
mselect=get(handles.navppselect,'Value');   
mselect=data.selmodes(mselect);             % what is the selection mode?
itrack=mselect>1;     % are we tracking?
data.imgf=(size(img,3)==3);                 % draw color or grayscale?
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);                   % is there editing-selection?
if(medt==8)
    % are we in proof-read mode, use alternative shading!!!
    shade=data.pref.shade(2);
else
    % use primary shading otherwise
    shade=data.pref.shade(1);
end

%=========================================================
% EXPERIMENTAL IMAGE REDUCTION MODE
sm=get(handles.axmain,'Position');
sm=sm([4,3])-sm([2,1]);
K=max(1,min(length(box{1})/sm(1),length(box{2})/sm(2)));
% EXPERIMENTAL IMAGE REDUCTION MODE
%==========================================================


% if flipped, no tracking
% SHADOW ADDON
if(data.shadow) itrack=0; end
if(iflip) itrack=0; end

% if necessary, update image data
if(iv1) utdgetimg(handles); end
if(iv2) utdgetidx(handles); end

% copy image to tmp3
tmp3=img;

% apply shades
if(itrack)
    if(data.imgf)
        tmp1=tmp3(box{:},:);
        ind=find(idx(box{:})); 
        sb=prod(size(idx(box{:})));
        ind=[ind(:);sb+ind(:);2*sb+ind(:)];
        tmp1(ind)=immultiply(tmp1(ind),shade);
        tmp3(box{:},:)=tmp1;                
        
        
%         for i=1:3
%             tmp=tmp3(box{:},i);
%             tmp(idx(box{:}))=immultiply(tmp(idx(box{:})),shade);
%             tmp3(box{:},i)=tmp;
%         end
    else
        if(~iedge)
            ind=find(idx(box{:}));
            tmp=tmp3(box{:});
            tmp(ind)=immultiply(tmp(ind),shade);
            tmp3(box{:})=tmp;
        else
            tmp3(1,1)=0; tmp3(1,2)=data.edgemax; 
        end
    end
end

% SHADOW ADDON
if(data.shadow & ~iflip)
    k=round(get(handles.navsldbrowse,'Value'));
    if(~iedge)
        if(strcmp(get(handles.uimshadowmbr,'Checked'),'off'))
            global cat
            tmp=tmp3(box{:});
            ind=find(cat{k}(box{:})>0 & ~idx(box{:}));
            tmp(ind)=immultiply(tmp(ind),data.pref.shade(2));
            ind=find(idx(box{:}));
            tmp(ind)=immultiply(tmp(ind),data.pref.shade(1));
            tmp3(box{:})=tmp;
        else
            global cat
            for i=1:3
                tmp=tmp3(box{:},i);
                ind=find(cat{k}(bow{:})>0 & ~idx(box{:}));
                tmp(ind)=immultiply(tmp(ind),data.pref.shade(2));
                ind=find(idx(box{:}));
                tmp(ind)=immultiply(tmp(ind),data.pref.shade(1));
                tmp3(box{:},i)=tmp;
            end
        end            
    else
        tmp3(1,1)=0; tmp3(1,2)=data.edgemax;
    end
end

% ACTUALLY draw
axes(handles.axmain);
%data.ih=imagesc(tmp3);
if(length(size(tmp3))<3 || size(tmp3,3)==1)
    if(~isa(tmp3,'uint8')) tmp3=im2uint8(tmp3); end
    data.ih=image(tmp3);
    clmp=gray(256); colormap(clmp);
else
    data.ih=image(tmp3);    
end
    
% Apply zoom
axis([data.zoom(3:4),data.zoom(1:2)]);
% Set colormap to gray if not color
colormap gray
% Make axis off so that can click in the image!!!
set(handles.axmain,'buttondownfcn',...
    'gui(''mainBtnDownFcn'',gcbo,[],guidata(gcbo))');
set(data.ih,'HitTest','off');

% alternatively::
%set(ih,'buttondownfcn',...
%        'get(get(gcbo,''parent''),''currentpoint'');');
%set(ih,'buttondownfcn',...
%    'guibrowser(''fmain_ButtonDownFcn'',gcbo,[],guidata(gcbo))');



% ######################################################################
function utdmove(k,handles)
% function performs shift in the stack to section # k
global al data proof

% this resets bwlabel used in quick-draw b/s new section
data.tmp{17}=[];

% implement move beyond stack boundaries
if(k==0 && isfield(proof,'prev') && exist(proof.prev))
    snext=proof.prev; dir=-1;
elseif(k>data.smax && isfield(proof,'next') && exist(proof.next))
    snext=proof.next; dir=1;
else
    snext=[]; dir=0; 
    % request current slice
    if(k>=1 & k<=length(al)) 
        ifget('alcat',k); 
    end
end

if(get(handles.navrbtcross,'value')) snext=[]; dir=0; end

if(~isempty(snext) && data.pref.ccdouble)
    data.cmove=data.cmove+1;
    if(data.cmove<3) return; end
end
data.cmove=0;

% implement move beyond stack boundaries
if(~isempty(snext) && exist(snext))
    global al cat debug
    set(handles.navtxtinfo,'String','Requesting adjacent substack...');
    set(handles.navtxtinfo,'fontweight','bold');
    pause(0.3);
    
    % identify global position
    a=k+proof.first-1;
    % identify global fov
    vzoom=data.zoom+proof.range([1,1,3,3])-1;
    % is there additional mapping?
    imap=get(handles.navrbtmap,'Value');
    % try to carry over the selection
    id=data.ids;
    if(dir==1) kadj=length(cat); else kadj=1; end    
    if(imap) 
        id=proof.pmap(id+1);
        [indx,indy]=find(ismember(proof.pmap(cat{kadj}(:,:)+1),id));
    else
        [indx,indy]=find(ismember(cat{kadj}(:,:),id));
    end    
    idsx=indx+proof.range(1)-1;
    idsy=indy+proof.range(3)-1;    
        
    
    % save current buffer
    global flMonitor
    if(isempty(flMonitor) || ~flMonitor)  
        % if monitoring, disable write-overs        
        s=sprintf('gui.buff%.3i.mat',proof.first);
        for k=1:length(cat) ind{k}=data.tmp{20}{k}; end
        
        % *** decided to disable saves of watershed-marks - expensive
%         if(isfield(debug,'ind1') && ~isempty(debug.ind1)) wcat=debug.ind1; 
%         else wcat=[]; end
        wcat=[];
    
         save(s,'-v6','ind','proof','wcat'); 
    else
        set(handles.navtxtinfo,'String',...
                        'monitoring - buffer not updated...');
        set(handles.navtxtinfo,'fontweight','bold');
        pause(0.3);
    end
    
    % load proof from adjacent substack
    tmp=proof;
    load(snext,'proof');
    
    % load cat & proof data
    if(~isfield(proof,'first'))
        % THIS IS FOR COMPATIBILITY
        set(handles.navtxtinfo,'String','Error loading adjacent substack');
        set(handles.navtxtinfo,'fontweight','bold');
        proof=tmp;
        return;
    else
        % THIS MANAGES UPDATES OF SUBSERIES FROM MEMORY BUFFERS
        global flUseMemBuffer        
        
        % loaded from buffer
        flLoadFromBuffer=0;

        % try to load buffer
        global buffer        
        if(~isempty(flUseMemBuffer) && flUseMemBuffer)
            debug.ind1=[];
            
            if(isempty(buffer)) buffer={}; end
            
            % try to reload the thing
            fstatus=0;
            while(~fstatus)
                iref=0;
                ffail=0;

                % see if we have copy in the buffer
                if(~isempty(buffer) && iscell(buffer))
                    for i=1:length(buffer)
                        if(~isempty(buffer{i}))
                            if(buffer{i}{1}==proof.first) iref=i; end
                        end
                    end
                end

                if(iref>0) % found - try to exchange
                    % buffer contains required data - swap
                    if(~isempty(flMonitor) && flMonitor) pause(0.25); end
                    pause(0.25);            
                    set(handles.navtxtinfo,'String','found memory buffer loading...');
                    set(handles.navtxtinfo,'fontweight','bold');            
                    pause(0.3);                                       
                    
                    buffer{iref}{1}=tmp.first;
                    step=1;
                    try
                        proof_old=tmp;
                        al_old=al;
                        cat_old=cat;
                        debug_old=debug;
                        ind_old=data.tmp{20};
                        step=2;

                        proof=buffer{iref}{2};
                        al=buffer{iref}{3};
                        cat=buffer{iref}{4};
                        debug=buffer{iref}{5};
                        data.tmp{20}=buffer{iref}{6};
                        
                        buffer{iref}{1}=proof_old.first;
                        buffer{iref}{2}=proof_old;
                        buffer{iref}{3}=al_old;
                        buffer{iref}{4}=cat_old;
                        buffer{iref}{5}=debug_old;
                        buffer{iref}{6}=ind_old;
                        
                        clear proof_old al_old cat_old debug_old ind_old
                        fstatus=1;
                        flLoadFromBuffer=1;
                    catch
                        if(step==2)
                            clear global al cat debug proof
                            clear al cat debug proof
                            
                            al=al_old;
                            cat=cat_old;
                            debug=debug_old;
                            data.tmp{20}=ind_old;

                            clear al_old cat_old debug_old ind_old
                        else
                            clear al_old cat_old debug_old ind_old
                        end
                        ffail=1;
                    end
                else % if not found, try to save anew
                    i=1;
                    while(i<=length(buffer) && ~isempty(buffer{i}))
                        i=i+1;
                    end

                    step=1;
                    try                        
                        % estimate reloading size
                        sz=whos('al','cat','debug','proof');                        
                        sz_sum=sum([sz.bytes]);
                        
                        l=0; sz_att=0;
                        for j=1:length(buffer) 
                            if(~isempty(buffer{j})) l=l+1; end; 
                            sz_att=max(sz_att,length(buffer{j}{6}));
                        end
                        sz_tmp20=0;
                        for j=1:length(data.tmp{20})
                            sz_tmp20=sz_tmp20+length(data.tmp{20}{j}{1});
                        end
                        sz_att=round(1.5*max(sz_att,sz_tmp20));
                        
                        sz_buf=whos('buffer');
                        if(l>0) sz_buf=sz_buf.bytes/l; else sz_buf=sz_sum; end
                        
                        % safe-buffer - need this to let sort below work
                        X=zeros(sz_att,1);
                        
                        % need this to prevent MATLAB going into upper 2G
                        %  will crash otherwise
                        Y={};
                        for j=1:20 Y{end+1}=zeros(1e6,1); end
                        
                        
                        % test-allocate memory
                        sz_est=sz(1).bytes*(1+sz_buf/sz_sum)/2*1.8/8/length(al);
                        al_old=cell(size(al));
                        for j=1:length(al) al_old{j}=zeros(round(sz_est),1); end
                        sz_est=sz(2).bytes*(1+sz_buf/sz_sum)/2*1.8/8/length(cat);
                        cat_old=cell(size(cat));
                        for j=1:length(cat) cat_old{j}=zeros(round(sz_est),1); end
                        sz_est=sz(3).bytes*(1+sz_buf/sz_sum)/2*1.8/8;
                        debug_old=zeros(round(sz_est),1);
                        sz_est=sz(4).bytes*(1+sz_buf/sz_sum)/2*1.8/8;
                        proof_old=zeros(round(sz_est),1);
                        
                        clear al_old cat_old debug_old proof_old Y
                        
                        buffer{i}={};
                        buffer{i}={tmp.first,tmp,al,...
                            cat,debug,data.tmp{20},0};
                        
                        step=2;

                        % reload stack from disk
                        set(handles.navtxtinfo,'String','loading...');
                        set(handles.navtxtinfo,'fontweight','bold');
                        pause(0.3);

                        
                        clear global al cat debug
                        clear al cat debug

                        global al cat debug
                        load(snext,'al','cat','debug');

                        % don't need this if loaded from buffer
                        s=sprintf('gui.buff%.3i.mat',proof.first);
                        if(exist(s)) load(s,'proof','ind'); wcat={}; 
                        else ind={}; wcat={}; end
                        clear X
                        fstatus=1;                        
                    catch
                        clear al_old cat_old debug_old proof_old Y
                        
                        if(step==2)
                            clear global al cat debug proof
                            clear al cat debug proof
                            proof=buffer{i}{2};
                            al=buffer{i}{3};
                            cat=buffer{i}{4};
                            debug=buffer{i}{5};
                            data.tmp{20}=buffer{i}{6};
                        end
                        ffail=1;
                        clear X
                    end
                end

                if(ffail) % try to free-up some space by purging buffer
                    ffail=0;
                    dref=-Inf; selfref=-1;
                    for i=1:length(buffer)
                        if(~isempty(buffer{i}))
                            if((abs(buffer{i}{1}-proof.first)+...
                                    abs(buffer{i}{1}-tmp.first))>dref)
                                dref=abs(buffer{i}{1}-proof.first)+...
                                    abs(buffer{i}{1}-tmp.first);
                                selfref=i;
                            end
                        end
                    end
                    
                    % there is nothing to purge
                    if(selfref==-1)
                        buffer={};
                        % reload stack from disk
                        clear global al cat debug
                        clear al cat debug

                        global al cat debug
                        load(snext,'al','cat','debug');

                        % don't need this if loaded from buffer
                        s=sprintf('gui.buff%.3i.mat',proof.first);
                        if(exist(s)) load(s);  else ind={}; wcat={}; end

                        fstatus=1;
                    else
                        buffer{selfref}={};
                        
                        X={}; cnt=0;
                        for i=1:length(buffer)
                            if(~isempty(buffer{i}))
                                cnt=cnt+1;
                                X{cnt}=buffer{i};
                            end
                        end
                        clear global buffer
                        clear buffer
                        global buffer
                        buffer=X;
                        clear X
                    end
                end
            end
        else
            % reload stack from disk
            clear global al cat debug
            clear al cat debug

            global al cat debug
            load(snext,'al','cat','debug');

            % don't need this if loaded from buffer
            s=sprintf('gui.buff%.3i.mat',proof.first);
            if(exist(s)) load(s);  else ind={}; wcat={}; end
        end
        
        if(~flLoadFromBuffer)
            data.tmp{20}=cell(size(cat));
            for k=1:length(ind)
                if(~isempty(ind{k}))
                    data.tmp{20}=ind{k};
                    cat{k}(ind{k}{1})=ind{k}{2};
                end                
            end
            if(~isempty(wcat)) debug.ind1=wcat; clear wcat; end                                             
        end
                
%         global fLoadInd2        
%         if(isempty(fLoadInd2) || fLoadInd2)
%             warning off MATLAB:load:variableNotFound
%             ind2=[];
%             try
%                 load(snext,'ind2');
%             catch
%                 fprintf('Failed loading ''ind2''\n');
%             end
%             debug.ind2=ind2; clear ind2
%             warning on MATLAB:load:variableNotFound
%         end        
        
        % reset color-mapping
        data.cmap=[];
    end

    % adjust zoom position
    vzoom([1,3])=max(proof.range([1,3]),vzoom([1,3]));
    vzoom([2,4])=min(proof.range([2,4]),vzoom([2,4]));
    vzoom=vzoom+1-proof.range([1,1,3,3]);     
    data.zoom=vzoom;
    
    k=a+1-proof.first;
    set(handles.navsldbrowse,'Value',k);
    
    data.nclst=0;
    data.logid=0;
    
    % reinitialize display & stuff
    gui_OpeningFcn(gcf,snext,handles);
    
    % tell user where we are
    if(k==1 || k==data.smax) s='*'; else s=''; end
    set(handles.navtxtslice,'String',...
                        sprintf('Section # %i%s',k+proof.first-1,s));
    
    % reset editing controls during cross-over
    set(handles.edtppcorrect,'Value',1);
    set(handles.edtcklink,'Value',0);
    set(handles.edtckdelete,'Value',0);
    set(handles.edtrbtdraw,'Value',0);
    set(handles.edtrbtquickdraw,'Value',0);
    set(handles.edtrbtdelete,'Value',0);
    data.tmp{4}=[];

    % update selection    
    indx=idsx-proof.range(1)+1;
    indy=idsy-proof.range(3)+1;
    
    if(k<=1)
        % I don't get what this is doing :(((
        set(handles.navtxtinfo,'String','Cannot maintain selection');
        set(handles.navtxtinfo,'fontweight','normal');
        data.ids=[];
    else
        kk=k-dir;
        indx=sub2ind(size(cat{kk}),indx,indy);
        indy=unique(cat{kk}(indx)); indy=setdiff(indy,0);
        
        if(imap) ids=proof.pmap(indy+1); else ids=indy; end
        ids=unique(ids);
        
        data.ids=ids;
        str='';
        for i=1:length(ids)  str=[str,sprintf('%i,',ids(i))]; end
        if(length(ids)>0) str=str(1:end-1); end
        set(handles.navedtselection,'String',str);                
    end    

    if(length(data.ids)==1) 
        % if have single obj in selection, make it also last selected
        data.nclst=data.ids;
        id=data.nclst;
    elseif(~isempty(data.edlogid) && ismember(data.edlogid(1),data.ids))
        % if have obj in selection identified in proof-note, make is also
        % last selected
        data.nclst=data.edlogid(1);
        id=data.nclst;
    else
        id=[];
    end

    if(~isempty(id))
        itype=proof.tmap(id+1);
        if(itype==1)
            itype=proof.tmap(proof.pmap==proof.pmap(id+1));
            itype=intersect(itype,[2,3,4,5]);
            if(length(itype)>1) itype=1; 
            elseif(length(itype)==0) itype=proof.tmap(id+1);  end
        end
        set(handles.edtpptype,'value',itype);
    else  set(handles.edtpptype,'value',1); end
        
    % update selection display
    utdgetidx(handles);
    utdshow(handles);
    
    set(handles.navtxtinfo,'String','done.');
    set(handles.navtxtinfo,'fontweight','bold');
    
    return
    % that's it    
end

% prevent move outside of the stack boundaries
k=max(1,min(data.smax,k));

% if we are at the top, disable further 'up' functions etc
set(handles.navbtup,'enable','on');
set(handles.navbtdown,'enable','on');
if(k==1 & (~isfield(proof,'prev') || ~exist(proof.prev))) 
    set(handles.navbtdown,'enable','off'); 
end
if(k==data.smax & (~isfield(proof,'next') || ~exist(proof.next))) 
    set(handles.navbtup,'enable','off'); 
end

% tell user where we are
if(k==1 || k==data.smax) s='*'; else s=''; end
set(handles.navtxtslice,'String',...
                        sprintf('Section # %i%s',k+proof.first-1,s));

% see if anything actually have changed and we need to redraw
kk=round(get(handles.navsldbrowse,'Value'));
if(kk==k) return; end
set(handles.navsldbrowse,'Value',k);

% redraw section
data.v1zoom=[];
data.v2zoom=[];
utdshow(handles);

% reset cross-substack move flag
data.cmove=0;

% update time stats
data.ctime(4)=data.ctime(4)+1;

% ####################################################################
function utdzoom(point,ffactor,handles)
% function changes zoom; center point 'point' and zoom factor 'ffactor'
% 0.5 - stays the same, 0.25 - zooms in and 1.0 - zooms out
global data

% get center point
c=min(data.shape,max([1,1],floor(point)));
% get current zoom
box=data.zoom;

% change zoom
for i=1:2
    w(i)=floor((box(2*i)-box(2*i-1))*ffactor);
    box(2*i-1)=c(i)-w(i);
    box(2*i)=c(i)+w(i);
end

% see to it that the stack boundaries are respected
w=2*w;
if(box(1)<1 | box(2)>data.shape(1))
    w(1)=min(w(1),data.shape(1)-1);
end
if(box(3)<1 | box(4)>data.shape(2))
    w(2)=min(w(2),data.shape(2)-1);
end
if(box(1)<1) box(1)=1; box(2)=1+w(1); end
if(box(1)>data.shape(1)) box(2)=data.shape(1); box(1)=box(2)-w(1); end

if(box(2)<1) box(1)=1; box(2)=1+w(1); end
if(box(2)>data.shape(1)) box(2)=data.shape(1); box(1)=box(2)-w(1); end

if(box(3)<1) box(3)=1; box(4)=1+w(2); end
if(box(3)>data.shape(2)) box(4)=data.shape(2); box(3)=box(4)-w(2); end

if(box(4)<1) box(3)=1; box(4)=1+w(2); end
if(box(4)>data.shape(2)) box(4)=data.shape(2); box(3)=box(4)-w(2); end
data.zoom=box;

% change zoom in the main axis, note inverted order
zoom=data.zoom;
v1zoom=data.v1zoom;
v2zoom=data.v2zoom;

%=============================================================
% EXPERIMENTAL ADDON FOR IMAGE REDUCTION
% adjustment due to image reduction
if(~isempty(v1zoom))
    x=v1zoom([2,4])-v1zoom([1 3]);
    y=zoom([2,4])-zoom([1 3]);
    if(max(2.5*y<x)) data.v1zoom=[]; v1zoom=[]; end
end

if(~isempty(v2zoom))
    x=v2zoom([2,4])-v2zoom([1 3]);
    y=zoom([2,4])-zoom([1 3]);
    if(max(2.5*y<x)) data.v2zoom=[]; v2zoom=[]; end
end
% EXPERIMENTAL ADDON FOR IMAGE REDUCTION
%=============================================================


if(isempty(v1zoom) | (zoom(1)<v1zoom(1)) | (zoom(2)>v1zoom(2)) |...
        (zoom(3)<v1zoom(3)) | (zoom(4)>v1zoom(4)) |...
        isempty(v2zoom) | (zoom(1)<v2zoom(1)) | (zoom(2)>v2zoom(2)) |...
        (zoom(3)<v2zoom(3)) | (zoom(4)>v2zoom(4)))
    % if out of processed labels, need to redraw
    utdshow(handles); 
else
    % otherwise may just change zoom
    axes(handles.axmain);
    axis([data.zoom(3:4),data.zoom(1:2)]);
end

% ####################################################################
function utdzoomin(handles,event)
% function realizes movein capacity, i.e. is zooming in on a
%   selected fragments to allow close viewing
global data debug tmp cat proof

if(nargin<2) event=0; end

if(~isfield(data,'pref')) w=100; else w=data.pref.w; end

idd=data.ids;
k=round(get(handles.navsldbrowse,'Value'));
ik=k-data.cpos+1;

imap=get(handles.navrbtmap,'Value');        % is there additional mapping?
if(imap)
    tmp=ismember(proof.pmap(cat{ik}+1),proof.pmap(idd+1));
else
    tmp=ismember(cat{ik},idd);
end

% if nothing found, use pmap anyway
if(isempty(find(tmp,1)))
    tmp=ismember(proof.pmap(cat{ik}+1),proof.pmap(idd+1));
end
w=round(1.5*w);

% extract bounding box for this fragment
stats=regionprops(uint8(tmp),'BoundingBox');

% if in proofreading mode tracking, check that object shifted by too-much
if(event & isempty(stats)) return; end

if(event)
    box=max(1,floor(stats(1).BoundingBox(end:-1:1)));
    for i=1:2
        vpos(2*i-1)=box(2+i);
        vpos(2*i)=box(2+i)+box(i);
    end
    dist=[vpos([1,3])-data.zoom([1,3]),data.zoom([2,4])-vpos([2,4])];
    dist=[min(dist([1,2])),min(dist([3,4]))];
    if(dist(1)>0.125*(data.zoom(2)-data.zoom(1)) &...
            dist(2)>0.125*(data.zoom(4)-data.zoom(3)))
        return;
    end
end

    
% unless there is nothing of this fragment in the slices...
if(~isempty(stats)) % change zoom    
    box=max(1,floor(stats(1).BoundingBox(end:-1:1)));
    for i=1:2
        data.zoom(2*i-1)=box(2+i)-w;
        data.zoom(2*i)=box(2+i)+box(i)+w;
    end
    data.zoom([1,2])=max(1,min(data.shape(1),data.zoom([1,2])));
    data.zoom([3,4])=max(1,min(data.shape(2),data.zoom([3,4])));
else
    set(handles.navtxtinfo,'String','No selection found');
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

    

% change zoom in the main axis, note inverted order
zoom=data.zoom;
v1zoom=data.v1zoom;
v2zoom=data.v2zoom;
if(isempty(v1zoom) | (zoom(1)<v1zoom(1)) | (zoom(2)>v1zoom(2)) |...
        (zoom(3)<v1zoom(3)) | (zoom(4)>v1zoom(4)) |...
        isempty(v2zoom) | (zoom(1)<v2zoom(1)) | (zoom(2)>v2zoom(2)) |...
        (zoom(3)<v2zoom(3)) | (zoom(4)>v2zoom(4)))
    utdshow(handles); 
else
    axes(handles.axmain);
    axis([data.zoom(3:4),data.zoom(1:2)]);
end


% #####################################################################
function mainBtnDownFcn(hObject,eventdata,handles)
% this function accepts clicks in the main axis
global data

ctype=get(gcf,'selectiontype');
cmode=data.pref.selmode;
fedt=get(handles.edtcklink,'Value') | ...
    get(handles.edtckdelete,'Value') | get(handles.edtrbtdraw,'Value') | ...
      get(handles.edtrbtdelete,'Value') | get(handles.edtrbtquickdraw,'Value') | ...
       get(handles.edtppcorrect,'Value')>1;
fedt=fedt & ~(~cmode & strcmp(ctype,'alt'));
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);                   % is there editing-selection?

 % click is edit click
if(fedt) edtmodify(handles); return; end

if(cmode)   % selection-biased clicking
    navmnselect_Callback(hObject,{0,ctype},handles);
    % if nontrivial editing mode, call edtselection
    if(ismember(medt,2:5)) edtselection(handles);  end
else        % context-menu biased clicking
    if(strcmp(ctype,'normal') | strcmp(ctype,'extend'))
        navmnselect_Callback(hObject,{1,ctype},handles);
        % if nontrivial editing mode, call edtselection
        if(ismember(medt,2:5)) edtselection(handles); end
    end
end


% ######################################################################
function navmnselect_Callback(hObject, eventdata, handles)
% handles on-click selection of objects in main axis
global proof data cat debug

% obtain point of click
point=get(handles.axmain,'currentpoint');
point=floor(point(1,end-1:-1:1));

if(min(point)<1 || max(point-data.shape)>0) return; end

% initialize:
% last selected id
data.nclst=0;

% obtain the slice # where selection is made
k=round(get(handles.navsldbrowse,'Value'));
ik=k-data.cpos+1;

% find out who've been clicked
id=double(cat{ik}(point(1),point(2)));
if(id==0 && data.shadow && ...
        strcmp(get(handles.uimshadowmbr,'checked'),'on') & ...
         isfield(debug,'ind2') && ~isempty(debug.ind2) && ...
         ~isempty(debug.ind2{ik}))
    id=double(debug.ind2{ik}(point(1),point(2)));
end

% is there additional mapping?
imap=get(handles.navrbtmap,'Value');        
if(imap) id=proof.pmap(id+1); end


% special behavior for accumulation of points
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);                   % is there editing-selection?
if(medt==10)
    global debug;
    if(strcmp(eventdata{2},'normal'))
        l=size(proof.points,2);
        proof.points(:,l+1)=[point(1),point(2),k,id]';
        set(handles.navtxtinfo,'String',...
            sprintf('select (%i,%i) %i',point(1),point(2),id));
        set(handles.navtxtinfo,'fontweight','normal');

        str=get(handles.dbgedtstats,'string');
        if(~isstruct(str) & strcmp(str,'none')) str={}; end
        l1=length(str);
        str{l1+1}=sprintf('select (%i,%i,%i) %i id# %i',point(1),point(2),k,id,l+1);
        set(handles.dbgedtstats,'string',str);

        % add red dot at the point of click
        range={point(1)-1:point(1)+1,point(2)-1:point(2)+1};
        range{1}=max(1,min(data.shape(1),range{1}));
        range{2}=max(1,min(data.shape(2),range{2}));
        debug.add1{k}(range{:})=255;

        utdgetimg(handles);
        utdshow(handles);
        return
    else
        global debug;
        set(handles.navtxtinfo,'String',...
            sprintf('remove (%i,%i)',point(1),point(2),id));
        set(handles.navtxtinfo,'fontweight','normal');
        i=(abs(proof.points(1,:)-point(1))<=1 & ...
            abs(proof.points(2,:)-point(2))<=1) & proof.points(3,:)==k;
        if(isempty(find(i)))
            set(handles.navtxtinfo,'String','Nothing to remove');
            set(handles.navtxtinfo,'fontweight','normal');
        else
            a=find(i);
            str=get(handles.dbgedtstats,'string');
            
            for j=a
                point=proof.points([1,2],i);
                range={point(1)-1:point(1)+1,point(2)-1:point(2)+1};
                range{1}=max(1,min(data.shape(1),range{1}));
                range{2}=max(1,min(data.shape(2),range{2}));
                debug.add1{k}(range{:})=0;
                
                l1=length(str);
                str{l1+1}=sprintf('removed (%i,%i,%i) X id# %i',point(1),point(2),k,j);                
            end
            
            proof.points=proof.points(:,~i);
            set(handles.dbgedtstats,'string',str);
            
            utdgetimg(handles);
            utdshow(handles);            
        end
        
        return
    end
end

% precaution for below selector
if(isempty(eventdata)) eventdata={1,'normal'}; end

% select action
altsel=data.pref.altsel;
switch(eventdata{1})
    case 0
        if(strcmp(eventdata{2},altsel{1}))   % regular select
            action=0;
        elseif(strcmp(eventdata{2},altsel{3}))  % box select
            action=2;
        elseif(strcmp(eventdata{2},altsel{2}))% regular deselect
            action=1;
        else
            action=-1;
        end
    case 1
        if(strcmp(eventdata{2},'normal'))   % regular select
            action=0;
        elseif(strcmp(eventdata{2},'extend'))  % regular deselect
            action=1;
        else
            action=-1;
        end
    otherwise
        action=-1;
end

% do what we have to do
% currently selected ids
str=get(handles.navedtselection,'string');
ids=data.ids(:)';
switch(action)
    case 0
        % left click - select
        % reset shift-click holder
        data.tmp{3}=[];
        % normal click, add to selection
        if(id>0)
            if(~ismember(id,ids))
                ids=[ids,id];
                if(isempty(str)) str=num2str(id);
                else str=[str,',',num2str(id)]; end
            end

            % process types:
            types=unique(proof.tmap(proof.pmap==proof.pmap(id+1)));
            types=intersect(types,2:5);
            if(length(types)>1)
                s='x';
            elseif(length(types)==1)
                s=get(handles.edtpptype,'String'); s=s{types}(1);
            else
                s='';
            end

            itype=proof.tmap(id+1);
            if(itype==1)
                itype=proof.tmap(proof.pmap==proof.pmap(id+1));
                itype=intersect(itype,[2,3,4,5]);
                if(length(itype)>1) itype=1; 
                elseif(length(itype)==0) itype=proof.tmap(id+1);  end
            end
            set(handles.edtpptype,'value',itype);

            data.nclst=id;
            data.ids=ids;
            set(handles.navedtselection,'String',str);
            set(handles.navtxtinfo,'string',sprintf('Selected %s%i',s,id));
            set(handles.navtxtinfo,'fontweight','normal');
        else
            % if nothing is selected - tell this to the operator
            set(handles.navtxtinfo,'string','None Selected');
            set(handles.navtxtinfo,'fontweight','normal');
        end
    case 1
        % right click - deselect
        % reset shift-select
        data.tmp{3}=[];
        
        if(ismember(id,ids))
            % remove id from selection
            ids=setdiff(ids,id);

            a=num2str(id); i=strfind(str,a);
            % cut id from selection string
            str=str([1:i-1,i+length(a)+1:end]);
            % remove possible ',' in the end
            if(~isempty(str))
                if(str(end)==',') str=str(1:end-1); end
                if(str(1)==',') str=str(2:end); end
            end

            % process types
            types=unique(proof.tmap(proof.pmap==proof.pmap(id+1)));
            types=intersect(types,2:5);
            if(length(types)>1)
                s='x';
            elseif(length(types)==1)
                s=get(handles.edtpptype,'String'); s=s{types}(1);
            else
                s='';
            end
            set(handles.edtpptype,'Value',1);            

            data.nclst=0;
            data.ids=ids;
            set(handles.navedtselection,'String',str);
            set(handles.navtxtinfo,'string',sprintf('Deselected %s%i',s,id));
            set(handles.navtxtinfo,'fontweight','normal');                        
        else
            set(handles.navtxtinfo,'string',sprintf('nothing to deselect'));
            set(handles.navtxtinfo,'fontweight','normal');
        end
    case 2
        % shift-click - box selection
        if(isempty(data.tmp{3}))
            data.tmp{3}=[point(1),point(2)];
            set(handles.navtxtinfo,'string',sprintf('marked UL'));
            set(handles.navtxtinfo,'fontweight','normal');            
        else
            ids=[];
            str=[];
            range={min(data.tmp{3}(1),point(1)):max(data.tmp{3}(1),point(1)),...
                min(data.tmp{3}(2),point(2)):max(data.tmp{3}(2),point(2))};
            id=unique(cat{ik}(range{:})); id=setdiff(id,union(0,ids));
            if(imap) id=unique(proof.pmap(id+1)); end
            id=id(:)';
            
            if(data.pref.selmajor) id=intersect(id,data.major); end                

            for i=id
                ids=[ids,i];
                if(isempty(str)) str=num2str(i);
                    else str=[str,',',num2str(i)]; end
            end
            data.tmp{3}=[];
            
            data.nclst=0;
            data.ids=ids;
            set(handles.navedtselection,'String',str);
            set(handles.edtpptype,'Value',1);              
            set(handles.navtxtinfo,'string',sprintf('selected'));
            set(handles.navtxtinfo,'fontweight','normal');            
        end
    otherwise
            set(handles.navtxtinfo,'string',sprintf('nothing to do'));
            set(handles.navtxtinfo,'fontweight','normal');
end

% call for debug-output if necessary, see 'setdebug' for more info
if(id>0 & data.pref.dbg)
    logid=dbgset(point(1),point(2),handles);
    if(~isempty(logid)) data.logid=logid;  else data.logid=0; end
else
    data.logid=0;
end

% redraw with new selection
utdgetidx(handles);
utdshow(handles);    


% ######################################################################
function navbtbackup_Callback(hObject, eventdata, handles)
% save backup button
savebackup(2,handles);

% ######################################################################
function navbtclear_Callback(hObject, eventdata, handles)
% handles clearing of selection on press of clear button
global proof data

% reset text-box holders of tracked clusters
% TODO: MAY NEED TO CHANGE AFTER ALL CONTROLS ARE IN
set(handles.navtxtinfo,'string','None Selected');
set(handles.navtxtinfo,'fontweight','normal');
set(handles.navedtselection,'string','');   
set(handles.edtpptype,'Value',1);
set(handles.edtppcorrect,'Value',1);
set(handles.edtcklink,'Value',0);
set(handles.edtckdelete,'Value',0);
set(handles.edtrbtdraw,'Value',0);
set(handles.edtrbtquickdraw,'Value',0);
set(handles.edtrbtdelete,'Value',0);

set(handles.dbgppevents,'Value',1,'String','none');
set(handles.dbgpplog,'Value',1,'String','none');
% set(handles.dbgckmajor,'Value',0);
set(handles.dbgedtstats,'String','none');

data.edconcl=0;     % conclusion made by operator for proof-note
data.edlogid=[];    % id of proof-note
data.edstats=[];    % miscellaneous statistics
data.ctime=zeros(1,7);% editing timing
tic;                % reset timer

% reset Editing controls state
set(handles.edtckcorrected,'Value',0);  % conclusion selector default
set(handles.edtcknotsure,'Value',0);    % conclusion selector default
set(handles.edtckother,'Value',0);      % conclusion selector default


% reset tracked clusters in data
data.nclst=0;
data.logid=0;
data.ids=[];
data.tmp{4}=[];
data.edlogid=[];

% redraw image with new selection
data.v2zoom=[];
utdgetidx(handles);
utdshow(handles);

% ######################################################################
function navbtdown_Callback(hObject, eventdata, handles)
% handles shift down in the stack on the click on 'down' button

% find out where we are
k=round(get(handles.navsldbrowse,'Value'))-1;
% actually move
utdmove(k,handles);

% ######################################################################
function navbtexit_Callback(hObject, eventdata, handles)
% handles click on exit button - correct exit from gui
global data proof

x=questdlg('Make sure to "save backup" before exiting!',...
    'FRIENDLY REMINDER...','Exit','Cancel','Exit');
if(~strcmp(x,'Exit')) return; end

% truncate proof.notes to only significant # of records
proof.notes=proof.notes(1:data.ednotes);

% if opened consoles - kill that
if(isfield(data,'fig1') && ~isempty(data.fig1)) delete(data.fig1); end
if(isfield(data,'fig2') && ~isempty(data.fig2)) delete(data.fig2); end

% deleting cross-substack buffers
if(data.pref.delbuf) delete gui*buff*; end


% SELF-DESTRUCT !!!
clear global data

% delete(handles.fmain);
closereq;

% ######################################################################
function navbtplay_Callback(hObject, eventdata, handles)
% supposed to handle playback in gui, now disabled
global data

if(~data.play)
    data.play=1;
    
    set(hObject,'String','stop');

    % find out where we are
    k=round(get(handles.navsldbrowse,'Value'));
    
    while(k<data.smax & data.play)
        k=min(data.smax,k+1);
        
        % actually move
        utdmove(k,handles);
        
        pause(0.33);
    end
else
    data.play=0;
    set(hObject,'String','Play');
end

data.play=0;
set(hObject,'String','Play');

% ######################################################################
function navbtreset_Callback(hObject, eventdata, handles)
% handles reseting of the zoom to full zoom on the click on 'reset' button
global proof data

% reset zoom for slide
data.zoom=[1,data.shape(1),1,data.shape(2)];
% apply new zoom
data.v1zoom=[];
data.v2zoom=[];

% change zoom in the main axis, note inverted order
zoom=data.zoom;
v1zoom=data.v1zoom;
v2zoom=data.v2zoom;
if(isempty(v1zoom) | (zoom(1)<v1zoom(1)) | (zoom(2)>v1zoom(2)) |...
        (zoom(3)<v1zoom(3)) | (zoom(4)>v1zoom(4)) |...
        isempty(v2zoom) | (zoom(1)<v2zoom(1)) | (zoom(2)>v2zoom(2)) |...
        (zoom(3)<v2zoom(3)) | (zoom(4)>v2zoom(4)))
    utdshow(handles); 
else
    axes(handles.axmain);
    axis([data.zoom(3:4),data.zoom(1:2)]);
end

% ######################################################################
function navbtup_Callback(hObject, eventdata, handles)
% handles shift up in the stack on the click on 'up' button
global proof data

% find out where we are
k=round(get(handles.navsldbrowse,'Value'))+1;
% actually move
utdmove(k,handles);

% ######################################################################
function navckoverlap_Callback(hObject, eventdata, handles)
% handles click on overlap radiobutton

utdgetimg(handles);
utdshow(handles);

% ######################################################################
function navckoverlay_Callback(hObject, eventdata, handles)
% overlay intialize function
global data

if(get(hObject,'Value'))
    % if clicked overlay, store first overlay participant
    mshow=get(handles.navppimgmode,'Value');    
    mshow=data.dismodes(mshow);                 % what is the display mode?    
    data.tmp{1}=mshow;
    if(data.tmp{1}<6)
        set(handles.navtxtinfo,'String','special overlay');
        set(handles.navtxtinfo,'fontweight','normal');
        data.tmp{1}=[];
        utdgetimg(handles);
        utdshow(handles);
    else
        set(handles.navtxtinfo,'String','overlay');   
        set(handles.navtxtinfo,'fontweight','normal');
    end
else
    % if unclicked, clear overlay things
    data.tmp{1}=[];
    utdgetimg(handles);
    utdshow(handles);
    set(handles.navtxtinfo,'String',{''});
    set(handles.navtxtinfo,'fontweight','normal');
end

% ######################################################################
function navedtgoto_Callback(hObject, eventdata, handles)
% implements quick goto
global proof data

% where we want to go
s=get(hObject,'String');
if(strcmp(s,'')) return; end
if(iscell(s)) s=s{1}; end
k=str2num(s)-proof.first+1;

if(k<1)
    set(handles.navtxtinfo,'String','manual move across substacks');
    set(handles.navtxtinfo,'fontweight','bold');
    k=1;
end

if(k>data.smax)
    set(handles.navtxtinfo,'String','manual move across substacks');
    set(handles.navtxtinfo,'fontweight','bold');
    k=data.smax;
end


% go there
if(~isempty(k)) utdmove(k,handles); end
% reset string value
set(hObject,'String',{''});
% TODO: MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);

% ######################################################################
function navedtselection_Callback(hObject, eventdata, handles)
% handles change of processes selection thorough edit-box
global proof data

% get the string from the text box
str=get(hObject,'String');
% identify separators
idd=findstr(str,',');
% if there is no ',' in the end, add by hands
if(~ismember(length(str),idd)) idd=[idd,length(str)+1]; end

% add first position to the list
idd=[1,idd];

% reset selection
ids=[];
% for each pair of positions, extract number and add it to the selection
for i=1:length(idd)-1
    str1=str(idd(i):idd(i+1)-1);
    ids=union(ids,str2num(str1));
end
% store selection in 'data', see description for data in the begining for
% more information
data.ids=ids;

% redraw with new selection
utdgetidx(handles);
utdshow(handles);
% TODO: MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);

% ######################################################################
function navmenu_Callback(hObject, eventdata, handles)


% ######################################################################
function navmnmove_Callback(hObject, eventdata, handles)
% handles move right-click function
global proof data

point=get(handles.axmain,'currentpoint');
point=point(1,end-1:-1:1);
utdzoom(point,0.5,handles);

% ######################################################################
function navmnmovein_Callback(hObject, eventdata, handles)
% handles movin right-click function
utdzoomin(handles);

% ######################################################################
function navmnzoomin_Callback(hObject, eventdata, handles)
% handles zoom in right-click function
point=get(handles.axmain,'currentpoint');
point=point(1,end-1:-1:1);
utdzoom(point,0.25,handles);

% ######################################################################
function navmnzoomout_Callback(hObject, eventdata, handles)
% handles zoom out right click function
point=get(handles.axmain,'currentpoint');
point=point(1,end-1:-1:1);
utdzoom(point,1.0,handles);

% ######################################################################
function navppimgmode_Callback(hObject, eventdata, handles)
% handles change of display mode
utdgetimg(handles);
utdshow(handles);
% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function navppselect_Callback(hObject, eventdata, handles)
% handles change of selection mode
utdgetimg(handles);
utdgetidx(handles);
utdshow(handles);
% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% --- Executes on button press in navrbtconsole.
function navrbtconsole_Callback(hObject, eventdata, handles)
% This function operates detached console.
global data

if(hObject==0)
    handles=guidata(data.fmain);
    hObject=handles.navrbtconsole;
    set(hObject,'Value',0);
end

if(data.pref.dbg) amax=155; bmax=0.75;
    else amax=125; bmax=0.98; end

if(get(hObject,'Value'))
    a=figure('Color',[0.7,0.7,0.7],'menubar','none',...
        'units','characters',...
        'position',[10 10 amax 12],'tag','fconsole','name','console',...
         'closerequestfcn','gui(''navrbtconsole_Callback'',0,[],[])');
    b=uicontrol('parent',a,'units','normalized',...
        'position',[0.01 0.05 bmax 0.90],'backgroundcolor',[1,1,1],...
        'enable','inactive','foregroundcolor',[0,0,0],...
        'max',2,'min',0,'string','',...
        'tag','dbgedtstats','style','edit','horizontalalignment','left');
    
    handles.dbgedtstatsold=handles.dbgedtstats;
    handles.fig1=a;
    handles.dbgedtstats=b;        

    if(data.pref.dbg)
        c=uicontrol('parent',a,'units','normalized',...
        'position',[0.76 0.8 0.24 0.1],'backgroundcolor',[1,1,1],...
        'enable','on','foregroundcolor',[0,0,0],...
        'max',1,'min',0,'string',{'none'},...
        'tag','dbgpplog','style','popupmenu','horizontalalignment','center');
        s='gcbo';
        s=['gui(''dbgpplog_Callback'',',s,',[],guidata(',s,'))'];
        set(c,'Callback',s);
        
         d=uicontrol('parent',a,'units','normalized',...
        'position',[0.76 0.6 0.24 0.1],'backgroundcolor',[1,1,1],...
        'enable','on','foregroundcolor',[0,0,0],...
        'max',1,'min',0,'string',{'none'},...
        'tag','dbgppevents','style','popupmenu','horizontalalignment','center');
        s='gcbo';
        s=['gui(''dbgppevents_Callback'',',s,',[],guidata(',s,'))'];
        set(d,'Callback',s);
    
        handles.dbgpplogold=handles.dbgpplog;
        handles.dbgppeventsold=handles.dbgppevents;
        
        handles.dbgpplog=c;
        handles.dbgppevents=d;
    end
    
    data.fig1=handles.fig1;
    data.fmain=handles.fmain;
    
    str=get(handles.dbgedtstatsold,'string');
    set(handles.dbgedtstats,'string',str);
    
    guidata(data.fig1,handles);
    if(~isempty(data.fig2)) guidata(data.fig2,handles); end
else
    if(isempty(data.fig1)) return; end
    handles.dbgedtstats=handles.dbgedtstatsold;
    if(isfield(handles,'dbgpplogold'))
        handles.dbgpplog=handles.dbgpplogold;
        handles.dbgppevents=handles.dbgppeventsold;
    end
    
    delete(handles.fig1);
    handles.fig1=[];
    data.fig1=[];
    guidata(data.fmain,handles);
    if(~isempty(data.fig2)) guidata(data.fig2,handles); end    
end

guidata(hObject,handles);

% ######################################################################
function navrbtflip_Callback(hObject, eventdata, handles)
% handles click on flip radiobutton
global data
data.flip=~data.flip;

% only need to redraw image with new selection
utdgetimg(handles);
utdshow(handles);

% ######################################################################
function navrbtmap_Callback(hObject, eventdata, handles)
% handles click on map radiobutton

% redraw image with map applied
utdgetimg(handles);
utdgetidx(handles);
utdshow(handles);

% --- Executes on button press in navrbtselstyle.
function navrbtselstyle_Callback(hObject, eventdata, handles)
% handles changes in clicking cycles
global data

% data.pref.selmode=get(hObject,'Value');
data.pref.selmode=~data.pref.selmode;
if(~isfield(data,'cntxt') || isempty(data.cntxt))
    data.cntxt=get(handles.axmain,'UIContextMenu');
end

if(data.pref.selmode)    
    set(handles.axmain,'UIContextMenu',[]);
else
%     handles=guidata(data.fmain);    
    set(handles.axmain,'UIContextMenu',data.cntxt);
end

% ######################################################################
function navsldbrowse_Callback(hObject, eventdata, handles)
% handles change of browser slider
global data

% get position and move
k=round(get(hObject,'Value'));
utdmove(k,handles);

% utdmove does not by itself redraw - redraw
data.v1zoom=[];
data.v2zoom=[];
utdshow(handles);

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function navsldmixer_Callback(hObject, eventdata, handles)
% handles change of mixer slider

% redraw image
utdgetimg(handles);
utdshow(handles);

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ======================================================================
%                           DEBUGING
% ######################################################################
function dbgedtstats_Callback(hObject, eventdata, handles)
% executes when changed content of dbgedtstats, nothing to do
% this control is empty


% ######################################################################
function dbgenable(mode,handles)
% function is responsible for enabling/disabling the bundle of
% debuging related controls

% TODO: CHANGE THESE AFTER ALL CONTROLS ARE IN

%set(handles.edtbtsave,'enable',mode);
set(handles.navppselect,'enable',mode);
%set(handles.edtppmode,'enable',mode);
%set(handles.navppimgmode,'enable',mode);
%set(handles.edtppcorrect,'enable',mode);
%set(handles.edtckmembrane,'enable',mode);
%set(handles.edtcklink,'enable',mode);
set(handles.navppselect,'enable',mode);
set(handles.edtbtcycle,'enable',mode);
set(handles.edtbtoverall,'enable',mode);
set(handles.edtcknotsure,'enable',mode);
%set(handles.edtckother,'enable',mode);
set(handles.edtckcorrected,'enable',mode);
set(handles.edtpptype,'enable',mode);
set(handles.edtbtprev,'enable',mode);
set(handles.edtbtnext,'enable',mode);
set(handles.edtppevents,'enable',mode);
%set(handles.edtppsize,'enable',mode);
set(handles.edtbtstroll,'enable',mode);

% ######################################################################
function idx=dbggetlogs(tag,nclst0,nclst,k,handles)
% returns logs relevant to fragment in slide k; input:
% tag is primary tag in slc;
% nclst0 is working label in cat0;
% nclst is global label in cat;
global proof data debug

if(~isfield(debug,'logs')) idx={}; return; end

% initialize output cell
idx=cell(3,1);

% find mergers pertainhig to selected cluster
% obtain log records associated with the cluster
idd=find(proof.pmap(1+debug.idd)==proof.pmap(1+nclst));
% select those that do not belong to lost
idx{2}=idd(~ismember(debug.ide(idd),[6,7]));

idd=find(proof.pmap(1+debug.idt)==proof.pmap(1+nclst));
idx{2}=union(idx{2},idd(ismember(debug.ide(idd),[6,7])));

% select events pertaining to this tag and also lookup for lost events 
% from the next slide by cat0 id
idd=find(debug.idk==k);

if(tag>0)
    % select all events refering this tag except [6,7] which refer cat0
    % and not slc! Will do [6,7] later
    idx{1}=idd(min([debug.logs(idd).tag]==tag,~ismember(debug.ide(idd),[6,7])));
    % add [6,7] events referenced by nclst0!!!
    idx{1}=union(idx{1},idd(min(ismember(debug.ide(idd),[6,7]),...
                                        [debug.logs(idd).tag]==nclst0)));
else
    idx{1}=idx{2}(debug.idk(idx{2})==k);
end

% adjust idx{2} to not include trivial events
idx{2}=idx{2}(~ismember(debug.ide(idx{2}),[2]));

% ######################################################################
function idd=dbggetmajor
% function returns the list of significant fragments in the proof
global cat data proof debug

% these are all processes that run through from first section to last
idd1=proof.pmap(unique(ifget('cat',1))+1);
idd2=proof.pmap(unique(ifget('cat',length(cat)))+1);
nidd=intersect(idd1,idd2);

switch(proof.v)
    % THIS SECTION IS FOR COMPATIBILITY:
    case {1,2,3}
        if(~isfield(debug,'logs')) idd=nidd; return; end
        ltmp=false([1,data.ednotes],'uint32');
        for i=1:data.ednotes ltmp(i)=~isempty(proof.notes(i).logid); end
        iddx=find(ltmp);
        
        idd1=[proof.notes(iddx).logid];        
        idd2=[proof.notes(iddx).concl];
        idd1=idd1(ismember(idd2,[2]));
        
        idd2=debug.ide(idd1);
        ltmp=ismember(idd2,[6,7]);
        nidd=union(nidd,debug.idt(idd1(ltmp)));
        nidd=union(nidd,debug.idd(idd1(~ltmp)));
        nidd=proof.pmap(nidd+1);        
    case {4,5,6}
        if(~isfield(debug,'logs')) idd=nidd; return; end
        ltmp=false([1,data.ednotes],'uint32');
        for i=1:data.ednotes ltmp(i)=~isempty(proof.notes(i).logid); end
        iddx=find(ltmp);
        
        idd1=[proof.notes(iddx).logid];        
        idd2=[proof.notes(iddx).concl];
        idd1=idd1(ismember(idd2,[1,3,5,7]));
        
        idd2=debug.ide(idd1);
        ltmp=ismember(idd2,[6,7]);
        nidd=union(nidd,debug.idt(idd1(ltmp)));
        nidd=union(nidd,debug.idd(idd1(~ltmp)));
        nidd=proof.pmap(nidd+1);        
    otherwise
        ltmp=false([1,data.ednotes]);
        for i=1:data.ednotes ltmp(i)=~isempty(proof.notes(i).tag); end
        iddx=find(ltmp);
        
        % all other processes that had been touched in proof-reading
        idd1=[proof.notes(iddx).tag];
        idd2=[proof.notes(iddx).concl];
        idd1=proof.pmap(idd1+1);
        %last record overrides everything
        stats=regionprops(idd1,'PixelIdxList','Area');
        ltmp=false(size(idd1));
        for i=1:length(idd1)
            ltmp(i)=ismember(idd2(stats(idd1(i)).PixelIdxList(end)),data.clist);
        end

        nidd=setdiff(union(nidd,idd1(ltmp)),idd1(~ltmp));
end
idd=intersect(nidd,proof.pmap(data.tmp{19}+1));

% ######################################################################
function logid=dbginfo(ntag,nclst0,nclst,k,handles)
% this prints log info specified by its ids
global proof debug cat data

if(~isfield(debug,'logs')) logid=0; return; end

% reset debuggin controls
set(handles.dbgpplog,'Value',1);
set(handles.dbgppevents,'Value',1);

% extract log records relevant to these ids
idx=dbggetlogs(ntag,nclst0,nclst,k,handles);

% fill in debug gui controls with log records
if(isempty(idx{1}))
    set(handles.dbgpplog,'String','none');
    logid=0;
else
    ids=zeros(size(idx{1}));
    str=cell(size(idx{1}));
    for i=1:length(idx{1})
        log=debug.logs(idx{1}(i));
        ids(i)=log.event;
        str{i}=['(',num2str(idx{1}(i)),') - ',num2str(log.event)];
    end
    set(handles.dbgpplog,'String',str);
    
    % smart selection of log record, depends on editing mode
    medt=get(handles.edtppmode,'Value');        
    medt=data.edtmodes(medt);           % is there editing-selection?
    switch(medt)
        case 2
            logid=find(ids==6);
        case 3
            logid=find(ids==3);
        case 4
            logid=find(ids==7);
        case 5
            logid=find(ids==9);
        otherwise
            logid=[];            
    end
    % position on the selected log record
    if(~isempty(logid))         
        set(handles.dbgpplog,'Value',logid(1));
        logid=idx{1}(logid(1));    
    else
        logid=idx{1}(1); 
    end
end

% fill out records in ppevents
if(isempty(idx{2}))
    set(handles.dbgppevents,'String','none');    
else    
    str=cell(size(idx{2}));
    for i=1:length(idx{2})
        log=debug.logs(idx{2}(i));
        str{i}=['(',num2str(idx{2}(i)),') - ',num2str(log.event),...
            ' @ ',num2str(log.slice)];
    end
    set(handles.dbgppevents,'String',str);
end

% print detailed information for current log record
if(logid>0) dbgstatout(logid,handles,ntag,nclst0,nclst); end

% ######################################################################
function dbgppevents_Callback(hObject, eventdata, handles)
% handles shift of focus on selection of log record from global log dropbox
global proof data debug

% get the position selected
id=get(hObject,'Value');
% get string representation of the position selected
str=get(hObject,'string');
if(iscellstr(str))  str=str{id}; end

% unless there is nothing selected...
if(~strcmp(str,'none'))
    % identify what is logid for this position
    pos=findstr(str,')');
    str=str(2:pos(1)-1);
    logid=str2num(str);
    data.logid=logid;
    
    % find out which slice corresponds to this record
    k=debug.logs(logid).slice;
    % if not there already, jump there
    if(k~=round(get(handles.navsldbrowse,'Value')))
        utdmove(k,handles);
    end
    % output stat info, see 'statout' in the utilities section for info
    dbgstatout(logid,handles);
    
    utdzoomin(handles);
end

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function dbgpplog_Callback(hObject, eventdata, handles)
% handles shift of focus by selecting different event from main log dropbox
global proof data

% get selected position and its string representation
id=get(hObject,'Value');
str=get(hObject,'string');
if(iscellstr(str)) str=str{id}; end

% unless nothing is in the box...
if(~strcmp(str,'none'))
    % find out what is logid for the record selected
    pos=findstr(str,')');
    str=str(2:pos(1)-1);
    logid=str2num(str);
    data.logid=logid;
    % output stat info, see 'statout' in the utilities section for info
    dbgstatout(logid,handles);    
end

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function logid=dbgset(x,y,handles)
% this finds logs for cluster selected by its point and outputs them
global proof debug cat data

% safeguard - if debug empty, return nothing
if(~isfield(debug,'logs')) logid=0; return; end

% extract segment id by its coordinates x,y,k
% if clicked membrane, bail out
k=round(get(handles.navsldbrowse,'Value'));
ik=k-data.cpos+1;
if(~isfield(debug,'slc') || isempty(debug.slc)) ntag=0; 
    else ntag=debug.slc{k}(x,y); end

% reference global label and working label from cat and cat0
nclst=cat{ik}(x,y);
if(isfield(debug,'cat0') && ~isempty(debug.cat0)) nclst0=debug.cat0{k}(x,y);
    else nclst0=nclst; end

% print debug info for selected set of ids
logid=dbginfo(ntag,nclst0,nclst,k,handles);

% ######################################################################
function dbgstatout(logid,handles,ntag,nclst0,nclst)
% function prints out detailed info for log record for #logid
global proof debug data

if(~isfield(debug,'logs')) return; end

% extract log record from debug.logs
log=debug.logs(logid);

% translate info from log record
tag=double(log.tag);
incom=double(log.incoming);             % incoming fragments
outgo=double(log.outgoing);             % outgoing fragment
ccidx=log.ccidx;                        % participagin fragments
flags=log.flags;                        % decision flags
ssums=log.ssums;                        % decision suppl. info
kk=length(ccidx);                       % # of participants
etype=log.event;                        % event type

% form output string
str=cell(1,kk+6);

% if no nclst cues had been provided, extract them from log record
if(nargin<5)
    switch(etype)
        case {6,7}
            ntag=0;
            nclst0=tag;
            nclst=double(debug.map(1+nclst0));
        otherwise
            ntag=tag;
            nclst0=outgo;
            nclst=double(debug.map(1+nclst0));
    end
end

% position onto the process by reference to its log record
data.nclst=nclst;
% output this process type /dendrite,axon,etc/
% set(handles.edtpptype,'Value',proof.tmap(1+nclst));
id=nclst;
itype=proof.tmap(id+1);
if(itype==1)
    itype=proof.tmap(proof.pmap==proof.pmap(id+1));
    itype=intersect(itype,[2,3,4,5]);
    if(length(itype)>1) itype=1; 
    elseif(length(itype)==0) itype=proof.tmap(id+1);  end
end
set(handles.edtpptype,'value',itype);

% pring-out detailed info
line=sprintf('[%i: %i, %i, %i(%i)], incoming:',...
    logid,ntag,nclst0,nclst,proof.pmap(1+nclst));
str{1}=line;
line=sprintf('%i,',incom);
line=[line(1:length(line)-1),...
    sprintf(' --> %i(%i)',outgo,proof.pmap(outgo+1))];
str{2}=line;

% print log record for each of the participating fragments
line=' id : flags : aux';
str{3}=line;

for i=1:kk
    % print cluster cat0 label
    line=sprintf('%5g:',ccidx(i));   
    
    % print decision flags
    if((size(flags,1)>0)&&(size(flags,2)>=i))
        line=[line,sprintf(' %3g |',flags(:,i))];
    end
    
    % print other suppl. parameters
    if((size(ssums,1)>0)&&(size(ssums,2)>=i))
        line=[line,sprintf(' %4.2f |',ssums(:,i))];
    end
    str{3+i}=line;
end
str{kk+4}='=================================';

% print note if any
if(~isempty(log.note))
    line=['note: ',sprintf(' %4.2f |',log.note)];
    str{5+kk}=line;
end

% print global statistics for this process
if(ismember(etype,[6,7])) nclst=debug.map(1+tag);
    else nclst=debug.map(1+outgo); end
if(isfield(debug,'sstats') && ~isempty(debug.sstats))
    line=['object 3D stats: ',sprintf(' %4.2f |',debug.sstats(nclst,:))];
end
str{kk+6}=line;
% write all this out
set(handles.dbgedtstats,'string',str);    


% =======================================================================
%                               EDITING 
% ######################################################################
function edtbtcycle_Callback(hObject, eventdata, handles)
% handles cycling through selected log record on the click in cycle button
global proof data

if(data.curfirst~=proof.first)
    set(handles.navtxtinfo,'string',...
        sprintf('you are in wrong substack, need %i',data.curfirst));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

% if in major mode, equivalent to 'prev'
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);                   % is there editing-selection?

% if not in editing mode, do nothing
if(isempty(data.tmp{2})) return; end

% if nothing had been selected so far, beging from record one...
flg=1;
if(isempty(data.fig2))
    fskip=get(handles.edtrbtthmode,'Value');
else
    fskip=strcmp(get(data.uimenu2(12),'Checked'),'on');
end

while(flg)
    % see if we need to keep scrolling
    str=get(handles.edtlbxevents,'String');
    iPos=get(handles.edtlbxevents,'Value');
    str=str{iPos};
        
    if(isempty(data.tmp{18})) % first step
        reload([],handles);
    else % otherwise get current position in ppInfo and shift it up
        iPos=get(handles.edtlbxevents,'Value');
        reload(iPos+1,handles);
    end
    
    iPos=get(handles.edtlbxevents,'Value');
    str=get(handles.edtlbxevents,'String');    
    str=str{iPos};
    flg=(fskip & strfind(str,'[1]') & (data.lpos<length(data.ind)));
end

% if we had flipped before, reset to nonflip state
data.flip=0;
set(handles.navppselect,'Value',find(data.selmodes==6));

% call edtppevents_Callback to finalize selection, 
%  see 'edtppevents_Callback' for info
edtppevents_Callback(handles.edtlbxevents,[],handles);

% ######################################################################
function edtbtnext_Callback(hObject, eventdata, handles)
% handles clicks on next button, see in the beginig for functionality
global proof data debug

if(data.curfirst~=proof.first)
    set(handles.navtxtinfo,'string',...
        sprintf('you are in wrong substack, need %i',data.curfirst));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

if(~isfield(data,'uislc') || isempty(data.uislc))
    return;
    kk=min(proof.nstack,round(get(handles.navsldbrowse,'Value'))+1);
else
    data.uislc=min(proof.nstack,data.uislc+1);
    kk=data.uislc;
    data.lpos=[];    
end
reload([],handles);
utdmove(kk,handles);


% ######################################################################
function edtbtprev_Callback(hObject, eventdata, handles)
% handles clicks on next button, see in the beginig for functionality
global proof data debug

if(data.curfirst~=proof.first)
    set(handles.navtxtinfo,'string',...
        sprintf('you are in wrong substack, need %i',data.curfirst));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

if(~isfield(data,'uislc') || isempty(data.uislc)) 
    return;
    kk=max(1,round(get(handles.navsldbrowse,'Value'))-1);
else
    data.uislc=max(1,data.uislc-1);
    kk=data.uislc;
    data.lpos=[];    
end
reload([],handles);
utdmove(kk,handles);

% #####################################################################
function edtbtsave_Callback(hObject, eventdata, handles)
% causes note to be saved if needed
edtsavenote(handles);

% ######################################################################
function edtbtstroll_Callback(hObject, eventdata, handles)
% handles clicks on 'stroll up' button
global data debug proof

% get what is the id in question
id=data.nclst;

if(~(id>0)) return; end

% find where it is lost
idd=find(proof.pmap==id)-1;
k=debug.objk(idd,2);
k=k(k<data.smax & k>0);
k=max(k);

% see what is the closest slice for this process and go there
utdmove(k,handles);

% do little trick with data.logid to zoomin on this fragment
utdzoomin(handles);

set(handles.navtxtinfo,'String','lost here...');
set(handles.navtxtinfo,'fontweight','normal');

% ######################################################################
function edtckcorrected_Callback(hObject, eventdata, handles)
% handles click on 'corrected' checkbox
global proof data debug

% save logid that is being currently active
% only do that once for a serries of clicks!
kk=round(get(handles.navsldbrowse,'Value'));
if(isempty(data.edlogid) || data.edlogid(1)==0) 
    if(data.nclst==0) set(hObject,'value',0); return; end
    data.edlogid=[data.nclst,kk];
    data.edfirst=proof.first;
    
    set(handles.navtxtinfo,'String',...
        sprintf('confirmed for %i @ %i',data.nclst,kk));
    set(handles.navtxtinfo,'fontweight','normal');
end

if(get(hObject,'value') && data.ctime(2)==0) data.ctime(2)=toc; end

% update editing timer
if(~get(handles.edtckcorrected,'value') && ...
        ~get(handles.edtcknotsure,'value') && ...
            ~get(handles.edtckother,'value'))
        data.ctime(2)=0;
        return;
end

% ######################################################################
function edtcklink_Callback(hObject, eventdata, handles)
% handles quick-link mode in & out
global data proof

if(get(hObject,'Value'))
    % currently selected clusters
    ids=unique(proof.pmap(data.ids+1));
    if(~isempty(data.edlogid) && data.edlogid(1)>0 && ...
            proof.first==data.edfirst) inclst=data.edlogid(1);
    else inclst=[]; end
    
    
    if(isempty(inclst))
        if(length(ids)==1)
            data.tmp{4}=ids;
            set(handles.navtxtinfo,'String',sprintf('quicklink to %i',ids));
            set(handles.navtxtinfo,'fontweight','normal');
            set(handles.edtppcorrect,'Value',2);
        else
            set(handles.navtxtinfo,'String','ambiguous: non single sel id!');
            set(handles.navtxtinfo,'fontweight','normal');
            set(hObject,'Value',0);
            data.tmp{4}=[];
        end
    else
        if((length(ids)>=1 & ...
                ismember(proof.pmap(inclst+1),ids)) | ...
                    isempty(ids))
            data.tmp{4}=inclst;
            set(handles.navtxtinfo,'String',sprintf('quicklink to %i',inclst));
            set(handles.navtxtinfo,'fontweight','normal');            
            set(handles.edtppcorrect,'Value',2);            
        else
            set(handles.navtxtinfo,'String',...
                'ambiguous: original proof id~=single sel id!');
            set(handles.navtxtinfo,'fontweight','normal');
            set(hObject,'Value',0);
            data.tmp{4}=[];
        end
    end
else
    set(handles.edtppcorrect,'Value',1);
    data.tmp{4}=[];
    set(handles.navtxtinfo,'String',sprintf('quicklink disabled'));
    set(handles.navtxtinfo,'fontweight','normal');
end

set(handles.edtckdelete,'Value',0);
set(handles.edtrbtdraw,'Value',0);
set(handles.edtrbtdelete,'Value',0);
set(handles.edtrbtquickdraw,'Value',0);

% ######################################################################
function edtckmembrane_Callback(hObject, eventdata, handles)
% handles quick-membrane mode
global data

set(handles.edtcklink,'Value',0); data.tmp{4}=[];
set(handles.edtrbtdraw,'Value',0);
set(handles.edtrbtdelete,'Value',0);
set(handles.edtrbtquickdraw,'Value',0);
set(handles.edtppcorrect,'Value',1);

% ######################################################################
function edtcknotsure_Callback(hObject, eventdata, handles)
% handles selection of 'notsure' checkbox
global data proof

% save logid that is being currently active
% only do that once for a serries of clicks!
kk=round(get(handles.navsldbrowse,'Value'));
if(isempty(data.edlogid) || data.edlogid(1)==0)
    if(data.nclst==0) set(hObject,'value',0); return; end
    data.edlogid=[data.nclst,kk];
    data.edfirst=proof.first;
    
    set(handles.navtxtinfo,'String',...
        sprintf('not sure for %i @ %i',data.nclst,kk));
    set(handles.navtxtinfo,'fontweight','normal');
else
    data.edlogid(2)=kk;
end

if(get(hObject,'value') && data.ctime(2)==0) data.ctime(2)=toc; end

% update editing timer
if(~get(handles.edtckcorrected,'value') && ...
        ~get(handles.edtcknotsure,'value') && ...
            ~get(handles.edtckother,'value'))
        data.ctime(2)=0;
        return;
end

% ######################################################################
function edtckother_Callback(hObject, eventdata, handles)
% handles selection of 'other' checkbox
global data proof

% save logid that is being currently active
% only do that once for a serries of clicks!
kk=round(get(handles.navsldbrowse,'Value'));
if(isempty(data.edlogid) || data.edlogid(1)==0)
    if(data.nclst==0) set(hObject,'value',0); return; end
    data.edlogid=[data.nclst,kk];
    data.edfirst=proof.first;
    
    set(handles.navtxtinfo,'String',...
        sprintf('other for %i @ %i',data.nclst,kk));
    set(handles.navtxtinfo,'fontweight','normal');
end

if(get(hObject,'value') && data.ctime(2)==0) data.ctime(2)=toc; end

% update editing timer
if(~get(handles.edtckcorrected,'value') && ...
        ~get(handles.edtcknotsure,'value') && ...
            ~get(handles.edtckother,'value'))
        data.ctime(2)=0;
        return;
end


% ######################################################################
function edtedtnotes_Callback(hObject, eventdata, handles)
% handles changes in notes edit-box, nothing for now

% ######################################################################
function edtenable(mode,handles)
% controls enabled/disabled bunch of controls for editing
%  select mode='on' or mode='off'

% TODO: AFTER ALL CONTROLS ARE IN
if(strcmp(mode,'on')) 
    dbgenable(mode,handles);
else
    set(handles.edtppcorrect,'enable',mode);
    set(handles.edtckmembrane,'enable',mode);
    set(handles.edtcklink,'enable',mode);
    set(handles.edtppcorrect,'enable',mode);
    set(handles.edtbtcycle,'enable',mode);
    set(handles.edtbtoverall,'enable',mode);
    set(handles.edtcknotsure,'enable',mode);
    set(handles.edtckcorrected,'enable',mode);
    set(handles.edtpptype,'enable',mode);
    set(handles.edtbtprev,'enable',mode);
    set(handles.edtbtnext,'enable',mode);
    set(handles.edtppevents,'enable',mode);
    set(handles.edtppsize,'enable',mode);
    
    set(handles.edtppmode,'Value',1);
    set(handles.edtppcorrect,'Value',1);
    set(handles.edtpptype,'Value',1);

    set(handles.edtedtnotes,'String','none');
    set(handles.edtcklink,'Value',0);
    set(handles.edtppevents,'String','none');
    set(handles.edtckmembrane,'Value',0);
    set(handles.edtppsize,'Value',1);
end

% ######################################################################
function [idd1,idd]=edtlfprep(k,handles)
% this function prepares the list of essential lost/found fragments that
% user will need to attend in the process
global proof debug data cat

global tmp

% find out what we need to do
imode=get(handles.edtppmode,'Value');
imode=data.edtmodes(imode);

% initialize log records pertaining to this slice
idd=find(debug.idk==k);
switch(imode)
    case 4 % find records for 'lost' events
        idd=idd(debug.ide(idd)==7);
        idd1=debug.idt(idd);
    case 5 % find records for 'found' events
        idd=idd(debug.ide(idd)==9);
        idd1=debug.idd(idd);
end

% bail out if none
if(isempty(idd)) return; end

% category II-III: globally important fragments, borrow from data.major
% identification will be made by global labels up to now
tmp=proof.pmap(1+idd1);
id1=ismember(tmp,proof.pmap(1+data.major));

% category IV: clusters must really be disappearing
% subtract clusters present in next/prev slice accordingly
switch(imode)
    case 4 
        if(k<data.smax) id1=id1 & ~ismember(tmp,proof.pmap(cat{k+1}+1)); end
    case 5 
        if(k>1) id1=id1 & ~ismember(tmp,proof.pmap(cat{k-1}+1)); end
end
idd=idd(id1);
idd1=idd1(id1);

id1=~ismember(idd1,proof.pmap(debug.esc));
idd=idd(id1);
idd1=idd1(id1);

% ######################################################################
function edtmodify(handles)
% function responds for clicks requesting slice modification when edtppmode
% and edtppcorrection are both nontrivial
global proof data cat debug

global tmp

% get the position of click
point=get(handles.axmain,'currentpoint');
point=max(1,floor(point(1,end-1:-1:1)));
% get the slice the click was made in
k=round(get(handles.navsldbrowse,'Value'));
ik=k-data.cpos+1;

% save backup if there had been already enough corrections
if(data.bkcount>=data.backup) savebackup(0,handles); end

% determine what exactly we do now:
if(get(handles.edtckdelete,'Value')) mode=1; end
if(get(handles.edtrbtdelete,'Value')) mode=2; end
if(get(handles.edtrbtdraw,'Value')) mode=4; end
if(get(handles.edtrbtquickdraw,'Value')) mode=3; end
if(get(handles.edtppcorrect,'Value')==2) mode=5; end
if(get(handles.edtppcorrect,'Value')==6) mode=6; end
if(get(handles.edtppcorrect,'Value')==7) mode=7; end
if(get(handles.edtppcorrect,'Value')==9) mode=8; end
if(get(handles.edtppcorrect,'Value')==5) mode=9; end
if(get(handles.edtppcorrect,'Value')==3) mode=10; end


switch(mode)
    case 1 % quick delete
        if(isempty(data.tmp{17}))
            data.tmp{17}=bwlabel(cat{ik}>0,data.pref.con);
        end
            
        tmp=data.tmp{17};
        col=tmp(point(1),point(2));
        ind=find(tmp==col);

        if(col==0)
            set(handles.navtxtinfo,'String','nothing to delete');
            set(handles.navtxtinfo,'fontweight','normal');
            return;
        end

        % delete
        tmp=cat{ik};        
        
        % undo buffer
        i=length(data.undo);
        data.undo(i+1).action=1;
        data.undo(i+1).slc=k;
        data.undo(i+1).data={tmp(ind),ind}; 
        data.undo(i+1).proof=proof.first;
        

        tmp(ind)=0;
        cat{ik}=tmp;        
        
        
        % remember in .tmp{20}
        ind=ind(:);
        if(isempty(data.tmp{20})) data.tmp{20}=cell(1,length(cat)); end
        
        if(isempty(data.tmp{20}{ik}))
            data.tmp{20}{ik}{1}=[]; 
            data.tmp{20}{ik}{2}=[];            
        end
        ltmp=ismember(data.tmp{20}{ik}{1},ind);
        ltmp1=ismember(ind,data.tmp{20}{ik}{1});
        data.tmp{20}{ik}{2}(ltmp)=cat{ik}(ind(ltmp1));
 
        ind=ind(~ltmp1);
        data.tmp{20}{ik}{1}=[data.tmp{20}{ik}{1}(:);ind];
        data.tmp{20}{ik}{2}=[data.tmp{20}{ik}{2}(:);cat{ik}(ind)];

        
        
        % set edtckcorrected to 1 to identify correction
        set(handles.edtckcorrected,'Value',1);
        edtckcorrected_Callback(handles.edtckcorrected,[],handles);
        utdgetidx(handles);        
        utdgetimg(handles);
        
        % update timing stats
        data.ctime(6)=data.ctime(6)+1;
    case 2 % regular delete
        range=cell(1,2);
        edrange=data.pref.edrange*2^(get(handles.edtppsize,'Value')-1);
        for i=1:2
            range{i}=max(1,point(i)-edrange):...
                min(data.shape(i),point(i)+edrange);
        end
        
        % undo buffer
        tmp=false(size(cat{ik}));
        tmp(range{:})=1;
        ind=find(tmp);
        tmp=cat{ik};        
        i=length(data.undo);
        data.undo(i+1).action=2;
        data.undo(i+1).slc=k;
        data.undo(i+1).data={tmp(ind),ind};
        data.undo(i+1).proof=proof.first;   
                  
        % delete
        cat{ik}(range{:})=0;
        
        % remember in .tmp{20}
        if(isempty(data.tmp{20})) data.tmp{20}=cell(1,length(cat)); end        
        
        if(isempty(data.tmp{20}{ik}))
            data.tmp{20}{ik}{1}=[]; 
            data.tmp{20}{ik}{2}=[];            
        end
        
        ind=ind(:);
        ltmp=ismember(data.tmp{20}{ik}{1},ind);
        ltmp1=ismember(ind,data.tmp{20}{ik}{1});
        data.tmp{20}{ik}{2}(ltmp)=cat{ik}(ind(ltmp1));
 
        ind=ind(~ltmp1);
        data.tmp{20}{ik}{1}=[data.tmp{20}{ik}{1}(:);ind];
        data.tmp{20}{ik}{2}=[data.tmp{20}{ik}{2}(:);cat{ik}(ind)];
        
        
        % force bwlabel if quick-edit again
        data.tmp{17}=[];

        % set edtckcorrected to 1 to identify correction
        set(handles.edtckcorrected,'Value',1);
        edtckcorrected_Callback(handles.edtckcorrected,[],handles);
        utdgetidx(handles);        
        utdgetimg(handles);
        
        % update timing stats
        data.ctime(7)=data.ctime(7)+1;        
    case 3 % quick draw
        if(~isempty(data.ids))
            ncol=data.ids(1);
        else
            set(handles.edtppcorrect,'Value',1);
            set(handles.navtxtinfo,'String','nothing to draw with');
            set(handles.navtxtinfo,'fontweight','bold');
            return
        end        
        
        if(isempty(data.tmp{17}))
            data.tmp{17}=bwlabel(cat{ik}>0,data.pref.con);
        end
            
        tmp=data.tmp{17};

        col=tmp(point(1),point(2));
        if(col==0)
            set(handles.navtxtinfo,'String','membrane selected');
            set(handles.navtxtinfo,'fontweight','normal');
            return;
        end        
        
        % undo buffer
        ind=find(tmp==col);
        tmp=cat{ik};
        i=length(data.undo);
        data.undo(i+1).action=3;
        data.undo(i+1).slc=k;
        data.undo(i+1).data={tmp(ind),ind};
        data.undo(i+1).proof=proof.first; 
        

        % recolor
        tmp(ind)=ncol;
        cat{ik}=tmp;        

        % remember in .tmp{20}
        if(isempty(data.tmp{20}{ik}))
            data.tmp{20}{ik}{1}=[]; 
            data.tmp{20}{ik}{2}=[];            
        end
        
        ind=ind(:);
        ltmp=ismember(data.tmp{20}{ik}{1},ind);
        ltmp1=ismember(ind,data.tmp{20}{ik}{1});
        data.tmp{20}{ik}{2}(ltmp)=cat{ik}(ind(ltmp1));
 
        ind=ind(~ltmp1);
        data.tmp{20}{ik}{1}=[data.tmp{20}{ik}{1}(:);ind];
        data.tmp{20}{ik}{2}=[data.tmp{20}{ik}{2}(:);cat{ik}(ind)];
        
        
        % set edtckcorrected to 1 to identify correction
        set(handles.edtckcorrected,'Value',1);
        edtckcorrected_Callback(handles.edtckcorrected,[],handles);        
        utdgetidx(handles);
        utdgetimg(handles);   
        
        % update timing stats
        data.ctime(6)=data.ctime(6)+1;        
    case 4 % regular draw
        if(~isempty(data.ids))
            ncol=data.ids(1);
        else
            set(handles.edtppcorrect,'Value',1);
            set(handles.navtxtinfo,'String','nothing to draw with');
            set(handles.navtxtinfo,'fontweight','bold');
            return
        end        
        
        range=cell(1,2);
        edrange=data.pref.edrange*2^(get(handles.edtppsize,'Value')-1);
        for i=1:2
            range{i}=max(1,point(i)-edrange):...
                min(data.shape(i),point(i)+edrange);
        end
        
        % undo buffer
        tmp=false(size(cat{ik}));
        tmp(range{:})=1;
        ind=find(tmp);
        tmp=cat{ik};
        i=length(data.undo);
        data.undo(i+1).action=2;
        data.undo(i+1).slc=k;
        data.undo(i+1).data={tmp(ind),ind};
        data.undo(i+1).proof=proof.first;
        

        % draw
        cat{ik}(range{:})=ncol;

        % remember in .tmp{20}
        if(isempty(data.tmp{20})) data.tmp{20}=cell(1,length(cat)); end
        
        if(isempty(data.tmp{20}{ik}))
            data.tmp{20}{ik}{1}=[]; 
            data.tmp{20}{ik}{2}=[];            
        end
        
        ind=ind(:);
        ltmp=ismember(data.tmp{20}{ik}{1},ind);
        ltmp1=ismember(ind,data.tmp{20}{ik}{1});
        data.tmp{20}{ik}{2}(ltmp)=cat{ik}(ind(ltmp1));
 
        ind=ind(~ltmp1);
        data.tmp{20}{ik}{1}=[data.tmp{20}{ik}{1}(:);ind];
        data.tmp{20}{ik}{2}=[data.tmp{20}{ik}{2}(:);cat{ik}(ind)];
        
        
        % force bwlabel if quick-edit again
        data.tmp{17}=[];        

        % set edtckcorrected to 1 to identify correction
        set(handles.edtckcorrected,'Value',1);
        edtckcorrected_Callback(handles.edtckcorrected,[],handles);
        utdgetidx(handles);        
        utdgetimg(handles);   
        
        % update timing stats
        data.ctime(7)=data.ctime(7)+1;        
    case 5 % this makes a new link
        % if tmp{4} is empty, this is our first click...
        if(isempty(data.tmp{4})) 
            % find out which cluster had been selected and store that in
            % tmp{4}
            data.tmp{4}=proof.pmap(1+double(cat{ik}(point(1),point(2))));
            
            % if membrane is selected, bail out!!!
            if((data.tmp{4}==0))
                set(handles.navtxtinfo,'String','MEMBRANE selected');
                set(handles.navtxtinfo,'fontweight','normal');
                set(handles.edtppcorrect,'Value',1);
                set(handles.edtcklink,'Value',0);
                data.tmp{4}=[];
                data.tmp{6}=[];
                return
            else
                set(handles.navtxtinfo,'String',...
                                ['Link ',num2str(data.tmp{4}),'--']);
                set(handles.navtxtinfo,'fontweight','normal');
                if(~ismember(data.tmp{4},proof.pmap(data.ids+1)))
                    data.ids=[data.ids,data.tmp{4}];
                    str=get(handles.navedtselection,'String');
                    str=[str,',',num2str(data.tmp{4})];
                    utdgetidx(handles);
                end
            end
        else  % if we are here, this is the second click
            % identify  the second cluster to be linked and store in tmp{6}
            data.tmp{6}=proof.pmap(1+double(cat{ik}(point(1),point(2))));
            
            % if either of these two is membrane, bail out
            if((data.tmp{4}==0) || (data.tmp{6}==0))
                set(handles.navtxtinfo,'String','MEMBRANE selected');
                set(handles.navtxtinfo,'fontweight','normal');
                data.tmp{6}=[];
                return
            else
                set(handles.navtxtinfo,'String',...
                    ['Link ',num2str(data.tmp{4}),' -- ',...
                     num2str(data.tmp{6})]);
                set(handles.navtxtinfo,'fontweight','normal'); 
            end
            
            
            % UPDATE TYPES            
            % process full types list
            idx1=unique(proof.tmap(proof.pmap==proof.pmap(1+data.tmp{4})));
            idx2=unique(proof.tmap(proof.pmap==proof.pmap(1+data.tmp{6})));            
            
            idx=union(idx1,idx2);
            idx=idx(ismember(idx,[2,3,4,5]));
            if(length(idx)>1)
                str=get(handles.navtxtinfo,'String');
                if(iscell(str)) str=str{1}; end
                str=[str,': type conflct'];
                set(handles.navtxtinfo,'String',str);
                set(handles.navtxtinfo,'fontweight','bold');
            else
                % update all n/a members on the first object's side
                idx=find(proof.pmap==proof.pmap(1+data.tmp{4}));
                idx=idx(ismember(idx1,1));
                i2=idx2(ismember(idx2,[2,3,4,5]));
                if(~isempty(i2))
                    % proof.tmap(idx)=i2;
                    set(handles.edtpptype,'Value',i2);
                end
                
                % update all n/a members on the second object's side
                idx=find(proof.pmap==proof.pmap(1+data.tmp{6}));
                idx=idx(ismember(idx2,1));
                i2=idx1(ismember(idx1,[2,3,4,5]));
                if(~isempty(i2))
                    % proof.tmap(idx)=i2;
                    set(handles.edtpptype,'Value',i2);
                end                                 
            end            
            
            % undo buffer
            idd=find(ismember(proof.pmap,[data.tmp{4},data.tmp{6}]));
            i=length(data.undo);
            data.undo(i+1).action=5;
            data.undo(i+1).slc=k;
            data.undo(i+1).data={proof.pmap(idd),idd};
            data.undo(i+1).proof=proof.first;            
                                  
            % check for previously attending the item
            ilast=ismember(proof.pmap(1+data.tmp{6}),...
                proof.pmap(1+[proof.notes(1:data.ednotes).tag]));
            ilast=max(double(ilast),data.smap(1+data.tmp{6}));            
            
            % make corrections to proof.pmap
            idd=find(proof.pmap==data.tmp{6});
            proof.pmap(idd)=proof.pmap(1+data.tmp{4});
                        
            % tell user where now it was last modified
            % ilast=data.smap(1+data.tmp{6});
            if(ilast>0)
                str=get(handles.navtxtinfo,'String');
                % if there is only one line, make it into cell string
                if(iscellstr(str)) str=str{1};  end

                str=[str,': last edt in '...
                    num2str(ilast)];
                set(handles.navtxtinfo,'string',str);
                set(handles.navtxtinfo,'fontweight','bold');
                pause(0.2);
            end
                                    
            
            % ALTERNATE THESE TO CHANGE BEHAVIOR OF quick link checkbox
            % if check box is on, allow for further quick-links
            if(~get(handles.edtcklink,'Value'))
                data.tmp{4}=[];                
                set(handles.edtppcorrect,'Value',1);
            end
            % reset other variables
            data.tmp{6}=[];        
            % set edtckcorrected to 1 to identify correction
            set(handles.edtckcorrected,'Value',1);
            % call edtckcorrected_Callback to do whatever else is necessary
            edtckcorrected_Callback(handles.edtckcorrected,[],handles);
            
            % reobtain shading map
            utdgetidx(handles);
            utdgetimg(handles);
            % tell GUI that significant objects may need to be recalculated
            data.vXupdt=1;
            
            % update timing stats
            data.ctime(5)=data.ctime(5)+1;
        end
    case 6 % restore clicked process to its original id#
        idds=cat{ik}(point(1),point(2));        
        idds=debug.map(idds+1);
        
        % make corrections to proof.pmap
        idd=find(debug.map==idds);
        
        % undo buffer
        i=length(data.undo);
        data.undo(i+1).action=6;
        data.undo(i+1).slc=k;
        data.undo(i+1).data={proof.pmap(idd),idd};
        data.undo(i+1).proof=proof.first;        
                                    
        proof.pmap(idd)=idds;
        
        % reset correction tools
        data.tmp{4}=[];
        set(handles.edtppcorrect,'Value',1);
        
        utdgetidx(handles);
        utdgetimg(handles);
        % tell GUI that significant object may need to be recalculated
        data.vXupdt=1;
        
        % update timing stats
        data.ctime(5)=data.ctime(5)+1;        
    case 7 % load last backup        
        savebackup(1,handles);
        
        utdgetidx(handles);
        utdgetimg(handles);
    case 8 % extend new ids
        i=length(proof.pmap);
        proof.pmap(end+1)=i;
        proof.tmap(end+1)=1;
        data.smap(end+1)=0;
        data.mmax=max(i,data.mmax);
        
        % make this id current selection
        data.ids=i;
        set(handles.navedtselection,'String',sprintf('%i',i));
        
        % make redraw image
        data.v1zoom=[];
        
        % remove any sort of selection
        global idx
        idx(:)=0;
        
        set(handles.navtxtinfo,'String',sprintf('pushed arrays to %i',i));
        set(handles.navtxtinfo,'fontweight','normal');
        data.tmp{4}=[];  
        set(handles.edtppcorrect,'Value',1);
    case 9
        edtundo(handles);
        set(handles.edtppcorrect,'value',1);
        data.tmp{4}=[];

        % update timing stats
        data.ctime(4)=data.ctime(4)+1;        
    case 10 % this groups objects
        if(isempty(data.ids)) return; end
        ids=data.ids;
        
        ids=unique(proof.pmap(ids+1));
        nid=min(ids);
        
        set(handles.navtxtinfo,'String',['Group into ',num2str(nid)]);
        set(handles.navtxtinfo,'fontweight','normal');
        

        % UPDATE TYPES
        % process full types list
        idx=unique(proof.tmap(ismember(proof.pmap,ids)));
        idx=idx(ismember(idx,[2,3,4,5]));
        if(length(idx)>1)
            str=get(handles.navtxtinfo,'String');
            if(iscell(str)) str=str{1}; end
            str=[str,': type conflcts'];
            set(handles.navtxtinfo,'String',str);
            set(handles.navtxtinfo,'fontweight','bold');
        elseif(length(idx)==1)
            % update all n/a members on the first object's side
            idx1=find(ismember(proof.pmap,ids));
            idx1=idx1(ismember(proof.tmap(idx1),1));
            if(~isempty(idx1))
                %proof.tmap(idx1)=idx;
                set(handles.edtpptype,'Value',idx);
            end
        end

        % undo buffer
        idd=find(ismember(proof.pmap,ids));
        i=length(data.undo);
        data.undo(i+1).action=6;
        data.undo(i+1).slc=k;
        data.undo(i+1).data={proof.pmap(idd),idd};
        data.undo(i+1).proof=proof.first;

        % make corrections to proof.pmap
        proof.pmap(idd)=nid;
                                    
            
        % set edtckcorrected to 1 to identify correction
        set(handles.edtckcorrected,'Value',1);
        % call edtckcorrected_Callback to do whatever else is necessary
        data.nclst=nid;
        edtckcorrected_Callback(handles.edtckcorrected,[],handles);

        % reobtain shading map
        utdgetidx(handles);
        utdgetimg(handles);
        % tell GUI that significant objects may need to be recalculated
        data.vXupdt=1;
        set(handles.edtppcorrect,'Value',1);
        
        % update timing stats
        data.ctime(5)=data.ctime(5)+1;        
end

% redraw to show corrections made
if(get(handles.edtrbtredraw,'Value')) utdshow(handles); end

% ######################################################################
function edtppcorrect_Callback(hObject, eventdata, handles)
% handles change of the editing modifications mode
global data

% clear link buffer
if(get(hObject,'Value')~=2) data.tmp{4}=[]; end

% see if we need to reset quick-links
if(get(handles.edtcklink,'Value'))
    data.tmp{4}=[];
    set(handles.edtcklink,'Value',0);
end
set(handles.edtckdelete,'Value',0);
set(handles.edtrbtdraw,'Value',0);
set(handles.edtrbtquickdraw,'Value',0);
set(handles.edtrbtdelete,'Value',0);

% only need to do something for watershed,load backup,group,UNDO and push id
if(ismember(get(hObject,'Value'),[3,5,7,9]))  edtmodify(handles); end

if(ismember(get(hObject,'Value'),10))
    set(hObject,'Value',1);        

    mkwatershed(handles);
end

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function edtppevents_Callback(hObject, eventdata, handles)
% handles selection of guided-inspection events from the list dropbox
global proof data debug proof

if(data.curfirst~=proof.first)
    set(handles.navtxtinfo,'string',...
        sprintf('you are in wrong substack, need %i',data.curfirst));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

% bail out if not in editing mode
if(isempty(data.tmp{2})) return; end

% see which log record we are requested to process
% need to scroll page?
str=get(handles.edtlbxevents,'String');
iPos=get(handles.edtlbxevents,'Value');
str=str{iPos};

flnext=0;
if(~isempty(findstr('next',str))) flnext=1;
elseif(~isempty(findstr('prev',str)))  flnext=-1; end

% if scroll page - need to reload - do not refocus
if(flnext) 
    reload([],handles); 
    
    if(flnext>0)
        set(handles.edtlbxevents,'Value',1);
    elseif(flnext<0)
        str=get(handles.edtlbxevents,'String');
        set(handles.edtlbxevents,'Value',length(str));
    end
    return;
end

% identify record for this entry
reload([],handles);
data.nclst=data.tmp{7}(data.ind(data.lpos),1);
dispslice=data.tmp{7}(data.ind(data.lpos),4);

% what slice needs to be shown
medt=get(handles.edtppmode,'Value');        
medt=data.edtmodes(medt);
if(medt==6) % in major mode slice id is different
    id=find(proof.pmap==proof.pmap(data.nclst+1))-1;
    kmin=min(data.objk(id,1));
    kmax=max(data.objk(id,2));
    if(kmin<=kmax && kmin>0) k=kmin;
     else set(handles.navtxtinfo,'String','Could not locate object'); end     
else
    k=data.tmp{5};
end

if(~isempty(data.uislc)) k=data.uislc; end

% if have prescribed location to show, do that
if(dispslice>0 & dispslice<=proof.nstack) k=dispslice; end

% if we are not there right now -- jump back
if(k~=round(get(handles.navsldbrowse,'Value')))
%     utdmove(k,handles);
    % facilitate move beyond stack boundaries
    if(k==1 & isfield(proof,'prev') && exist(proof.prev))
        snext=proof.prev; dir=-1;
    elseif(k==data.smax & isfield(proof,'next') && exist(proof.next))
        snext=proof.next; dir=1;
    else
        snext=[]; dir=0;
    end
    
    set(handles.navsldbrowse,'Value',k);
    if(k==1 || k==data.smax) s='*'; else s=''; end
    set(handles.navtxtslice,'String',...
        sprintf('Section # %i%s',k+proof.first-1,s));
    utdgetimg(handles);

    set(handles.navbtup,'enable','on');
    set(handles.navbtdown,'enable','on');
    if(isempty(snext))
        if(k==1) set(handles.navbtdown,'enable','off'); end
        if(k==data.smax) set(handles.navbtup,'enable','off'); end
    end
end

% this is a trick for edtselection to work properly 
set(handles.navppselect,'Value',find(data.selmodes==6));

% call 'editselection' for futher processing, see 'editselection' for info
edtselection(handles);

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function edtppmode_Callback(hObject, eventdata, handles)
% handles editing information change requests
global proof data debug cat tmp

% retrieve editing mode requested
iEdit=get(hObject,'Value');
iEdit=data.edtmodes(iEdit);

% retrieve current slide position
k=round(get(handles.navsldbrowse,'Value'));

% check we have something to work with
if(~ismember(iEdit,[1,6,8,9,10]) & ~isfield(debug,'logs'))
    set(handles.navtxtinfo,'String','no debug record...');
    set(handles.navtxtinfo,'fontweight','normal');
    set(hObject,'Value',1);
    return;
end

% reset all important variables if editing
data.tmp{2}=[];     % current list of to-do records
data.tmp{4}=[];     % current link member for quick-link
data.tmp{7}=[];     % state of proofing list
data.tmp{11}=[];    % list of not sure records
data.tmp{17}=[];    % bwlabel used for quick-draw
data.tmp{18}=[];    % last select obj id from to-do list

if(~ismember(iEdit,[1,7]))
    %data.edfirst=proof.first;
    data.edfirst=[];
    data.modefirst=proof.first;

    data.edconcl=0;     % conclusion made by operator for proof-note
    data.edlogid=[];    % id of proof-note
    data.edstats=[];    % miscellaneous statistics

    data.logid=0;       % currently select log record
    data.cslc=k;        % currently processed slice
    data.nclst=0;       % currently selected object    
    
    data.ids=[];        % set of selected objects    

    % reset Editing controls state
    set(handles.edtckcorrected,'Value',0);  % conclusion selector default
    set(handles.edtcknotsure,'Value',0);    % conclusion selector default
    set(handles.edtckother,'Value',0);      % conclusion selector default
    set(handles.edtpptype,'Value',1);       % process type to default

    set(handles.edtedtnotes,'String','');       % notes text box default
    set(handles.edtbtcycle,'enable','on');      % cycle button on
    set(handles.edtppcorrect,'Value',1);        % correct dropbox default
    set(handles.edtppsize,'Value',2);           % pen-size default
    set(handles.edtlbxevents,'Value',1);        % position in list-dropbox
    set(handles.navppselect,'Value',2);         % tracking
    set(handles.navrbtmap,'Value',1);           % mapping

    set(handles.edtcklink,'Value',0);           % quick link to default
    set(handles.edtckdelete,'Value',0);         % quick membrane to default
    set(handles.edtrbtdraw,'Value',0);          % quick draw to default        
    set(handles.edtrbtdelete,'Value',0);        % quick draw to default

    set(handles.navedtselection,'String','');      % selection string
    set(handles.navtxtinfo,'String','');        % info string
    set(handles.navtxtinfo,'fontweight','normal');    
    
    set(handles.edtbtnext,'Enable','Off');       % enable next button
    set(handles.edtbtprev,'Enable','Off');      % enable prev button
    
    data.bkcount=data.backup;                   % force backup save at 
                                                % first editing                                                
else    % if not editing, disable cycle and overall
    data.edfirst=[];    
    data.modefirst=[];

    % reset Editing controls state
    set(handles.dbgedtstats,'String',{''}); % reset messages console    
    set(handles.edtckcorrected,'Value',0);  % conclusion selector default
    set(handles.edtcknotsure,'Value',0);    % conclusion selector default
    set(handles.edtckother,'Value',0);      % conclusion selector default

    set(handles.edtppcorrect,'Value',1);        % correct dropbox default
    set(handles.edtcklink,'Value',0);           % quick link to default
    set(handles.edtckdelete,'Value',0);         % quick membrane to default
    set(handles.edtrbtdraw,'Value',0);          % quick draw to default        
    set(handles.edtrbtdelete,'Value',0);        % quick draw to default
    set(handles.edtrbtthmode,'value',0);    % set default through mode    

    
    set(handles.edtbtnext,'Enable','Off');       % enable next button
    set(handles.edtbtprev,'Enable','Off');      % enable prev button
    
    
    set(handles.edtbtcycle,'enable','off');
    set(handles.edtlbxevents,'Value',1,'String','none');
    
    set(handles.edtbtnext,'Enable','Off');       % enable next button
    set(handles.edtbtprev,'Enable','Off');      % enable prev button    
end

switch(iEdit)
    case 2 % if asked to edit splits, select split events for this slice
        idd=find(debug.idk==k);
        idd=idd(debug.ide(idd)==6);
        idd=debug.idt(idd);
    case 3 % if asked to edit mergers, select merger events for this slice
        idd=find(debug.idk==k);
        idd=idd(debug.ide(idd)==3);
        idd=debug.map(debug.idd(idd)+1);
    case 4 % if asked to edit losts, call for edtlfprep, see 'edtlfprep'
        idd=edtlfprep(k,handles);
    case 5 % if asked to edit founds, call for edtlfprep, see 'edtlfprep'
        idd=edtlfprep(k,handles);
    case 6 % mode for all-major processes
        idd=data.major;
    case 9 % extract 'not sure' notes
        set(handles.navtxtinfo,'String','processing notes...');
        set(handles.navtxtinfo,'fontweight','normal');
        pause(0.05);
        
        idd=[proof.notes(1:data.ednotes).concl];
        idd=find(ismember(idd,data.nslist));
        idd=idd([proof.notes(idd).slice]==k);
        
        data.tmp{11}=idd;        
        idd=[proof.notes(idd).tag];        

        set(handles.navtxtinfo,'String','');
        set(handles.navtxtinfo,'fontweight','normal');
        set(handles.dbgedtstats,'string','');
    case 10
        if(data.allowgather)
            proof.points=[];
            if(~isfield(debug,'add1') || isempty(debug.add1))
                debug.add1=cell(size(cat));
                for l=1:length(cat) 
                    debug.add1{l}=zeros(size(cat{l}),'uint8');
                end
            end
            idd=[];
            set(handles.navrbtmap,'Value',0);
            set(handles.dbgedtstats,'String','');
            set(handles.edtbtcycle,'enable','off');
            % initialize overlay
            set(handles.navckoverlay,'Value',1);
            data.tmp{1}=6;
            set(handles.navppimgmode,'Value',find(data.dismodes==13));
            set(handles.navppselect,'Value',1);
            utdgetimg(handles);
            utdgetidx(handles);
            utdshow(handles);
        else % if gather is forbiden, refuse to do anything
            set(handles.edtppmode,'Value',1);
        end

        return;
    otherwise % otherwise reset to inactive state
        data.tmp{2}=[];
        data.tmp{7}=[];
        data.tmp{11}=[];
        data.cslc=0;
        idd=[];        
end

% % if enabled through-mode -- reduce by processed clusters
% if(ismember(iEdit,[2,3,4,5,6]) & get(handles.edtrbtthmode,'Value'))
%     % find all finalized objects
%     id1=[proof.notes(1:data.ednotes).concl];
%     id1=find(ismember(id1,data.clist));
%     id1=proof.pmap([proof.notes(id1).tag]+1);
%     % only allow discarding of specific object types
%     id1=id1(~ismember(proof.tmap(id1+1),data.pref.objdskt));
%     
%     % find all elements of finilized objects
%     ltmp=ismember(proof.pmap(idd+1),id1);
% 
%     idd=idd(~ltmp);
%     if(~isempty(idm)) idm=idm(~ltmp); end
% end    

if(~ismember(iEdit,[1,7]))
    set(handles.edtrbtthmode,'value',1);    % set default through mode
else
    set(handles.edtrbtthmode,'value',0);    % set default through mode
end

% store selected objs ids
if(iEdit==6)
    id1=union([proof.notes(1:data.ednotes).tag],proof.pmap(idd+1));
    id1=setdiff(id1,idd); idd=idd(:); id1=id1(:);
    
    idx=[idd;id1];

    data.tmp{7}=zeros(length(idx),4);
    data.tmp{7}(:,1)=idx;
    % locations of data.major in extended unique list
    idx=1:length(idd);
    data.tmp{7}(idx,2)=data.morder(:);
    data.tmp{7}(idx,3)=data.vorder(:);
    data.tmp{7}(idx,4)=data.zmajor(:);
    % data.tmp{7} :: obj ID  |  obj ordering  |   obj volume | z-loc ::
    
    % this holds filter/marking options
    data.tmp{8}=false(size(data.tmp{7},1),4);
    % data.tmp{8} :: obj filter | obj NS | obj checked | obj shown ::

    % reassign the list of objects
    idd=data.tmp{7}(:,1);
    
    % reset initial position
    data.lpos=[];
    data.ind=[];
else
    idd=idd(:);
    % force reloading of the proofing list
    data.tmp{7}=[];
end

data.tmp{2}=idd;
% currently active slice
data.tmp{5}=k;
% corrected records
data.tmp{9}=zeros(size(idd));
% set indicator for which substack is being edited
data.curfirst=proof.first;
% set update major list by default
data.vXupdt=1;

% intialize display of selected log records in edtppevents dropbox
set(handles.edtlbxevents,'UserData',[-1 -1]);
data.lpos=[];
reload([],handles);
set(handles.edtlbxevents,'Value',1);

% redraw with new selection
utdgetimg(handles);
utdgetidx(handles);
utdshow(handles);

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function edtppsize_Callback(hObject, eventdata, handles)
% handles changes of pen size, nothing for now

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function edtpptype_Callback(hObject, eventdata, handles)
% handles selection of process type in the type dropbox
global proof data

% if a process being selected
if(data.nclst>0)    
    itype=get(handles.edtpptype,'Value');
    itype_old=proof.tmap(1+data.nclst);    
    if(ismember(itype,1:5)) 
        % reset all other primary type carriers
        idx=find(proof.pmap==proof.pmap(1+data.nclst));
        idx=idx(ismember(proof.tmap(idx),2:5));
        proof.tmap(idx)=1; 
    else
        % if changed id is primary type carrier & chosen secondary type
        if(ismember(itype_old,2:5))
            idx=find(proof.pmap==proof.pmap(1+data.nclst));
            % if no other primary type carriers found
            if(isempty(find(ismember(proof.tmap(idx),2:5),1)))
                % find an empty type id, if any
                idx=idx(proof.tmap(idx)==1);
                % assign it as primary type carrier
                if(~isempty(idx)) proof.tmap(idx(1))=itype_old; end
            end
        end
    end
    % reassign type
    proof.tmap(1+data.nclst)=itype;

    
%     ind=find(proof.pmap==proof.pmap(1+data.nclst));
%     idx=proof.tmap(ind);
%     % if assigning a major type - only affect major types
%     if(ismember(itype,[1,2,3,4,5]))
%         proof.tmap(ind(ismember(idx,[1,2,3,4,5])))=itype;
%         proof.tmap(1+data.nclst)=itype;
%     end
%     % if assigning a submajor type - only affect n/a's
%     if(ismember(itype,[1,6,7,8,9,10,11]))
%         proof.tmap(ind(ismember(idx,1)))=itype;
%         proof.tmap(1+data.nclst)=itype;
%     end
        
    % update type string
    id=data.nclst;
    types=unique(proof.tmap(proof.pmap==proof.pmap(data.nclst+1)));
    types=intersect(types,2:5);
    if(length(types)>1)
        s='x';
    elseif(length(types)==1)
        s=get(handles.edtpptype,'String'); s=s{types}(1);
    else
        s='';
    end
    set(handles.navtxtinfo,'string',sprintf('Selected %s%i',s,id));
    set(handles.navtxtinfo,'fontweight','normal');
else    % otherwise reset
    set(handles.edtpptype,'Value',1);
end

% MAKE SURE THE CONTROL LOOSES FOCUS
uicontrol(handles.navbtreset);


% ######################################################################
function edtrbtdelete_Callback(hObject, eventdata, handles)
% handles quick-membrane mode
global data

set(handles.edtcklink,'Value',0); data.tmp{4}=[];
set(handles.edtckdelete,'Value',0);
set(handles.edtrbtdraw,'Value',0);
set(handles.edtrbtquickdraw,'Value',0);
set(handles.edtppcorrect,'Value',1);

function edtrbtdraw_Callback(hObject, eventdata, handles)
% hanldes quick fragment drawing link
global data

set(handles.edtcklink,'Value',0); data.tmp{4}=[];
set(handles.edtckdelete,'Value',0);
set(handles.edtrbtdelete,'Value',0);
set(handles.edtrbtquickdraw,'Value',0);
set(handles.edtppcorrect,'Value',1);

% --- Executes on button press in edtrbtproofing.
function edtrbtproofing_Callback(hObject, eventdata, handles)
% this controls proofing window

% This function operates detached console.
global data

if(hObject==0)
    handles=guidata(data.fmain);
    hObject=handles.edtrbtproofing;
    set(hObject,'Value',0);
end

if(get(hObject,'Value'))
    a=figure('Color',[0.7,0.7,0.7],'menubar','none',...
        'units','characters',...
        'position',[10 10 60 40],'tag','fproof','name','proofing',...
         'closerequestfcn','gui(''edtrbtproofing_Callback'',0,[],[])');
    b=uicontrol('parent',a,'units','normalized',...
        'position',[0.02 0.02 0.95 0.96],'backgroundcolor',[1,1,1],...
        'enable','on','foregroundcolor',[0,0,0],...
        'max',1,'min',0,'string','',...
        'tag','edtlbxevents','style','listbox',...
        'horizontalalignment','left','FontName','FixedWidth');
    set(b,'Callback','gui(''edtppevents_Callback'',gcbo,[],guidata(gcbo))');
    
    d=[];
    c1=uimenu('Label','Filter','Tag','uifilter');
    d=uimenu(c1,'Label','Major','Tag','uifiltermajor','Checked','on',...
        'Callback','gui(''uifilter_Callback'',gcbo,1,guidata(gcbo))');
    d=[d,uimenu(c1,'Label','Checked','Tag','uifilterchecked',...
        'Callback','gui(''uifilter_Callback'',gcbo,2,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','Not Checked','Tag','uifilterunchecked',...
        'Callback','gui(''uifilter_Callback'',gcbo,3,guidata(gcbo))')];        
    d=[d,uimenu(c1,'Label','Has Type','Tag','uifilterchecked',...
        'Callback','gui(''uifilter_Callback'',gcbo,4,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','No Type','Tag','uifilterunchecked',...
        'Callback','gui(''uifilter_Callback'',gcbo,5,guidata(gcbo))')];            
    d=[d,uimenu(c1,'Label','Not Sure','Tag','uifilterNS',...
        'Callback','gui(''uifilter_Callback'',gcbo,6,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','Other','Tag','uifilterOther',...
        'Callback','gui(''uifilter_Callback'',gcbo,7,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','By Slice','Tag','uifiltercurrent',...
        'Callback','gui(''uifilter_Callback'',gcbo,8,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','By Type','Tag','uifiltertype',...
        'Callback','gui(''uifilter_Callback'',gcbo,9,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','AUTO','Tag','uifilterauto',...
        'Callback','gui(''uifilter_Callback'',gcbo,10,guidata(gcbo))')];
    
        
    c2=uimenu('Label','Sorting','Tag','uisort');
    d=[d,uimenu(c2,'Label','Sort Volume','Tag','uisortplain','Checked','on',...
        'Callback','gui(''uifilter_Callback'',gcbo,11,guidata(gcbo))')];
    
    c3=uimenu('Label','Other','Tag','uiother');
    d=[d,uimenu(c3,'Label','Skip Done','Tag','uiskip','Checked','on',...
        'Callback','gui(''uifilter_Callback'',gcbo,12,guidata(gcbo))')];
    d=[d,uimenu(c3,'Label','P-Map','Tag','uipmap','Checked','on',...
        'Callback','gui(''uifilter_Callback'',gcbo,13,guidata(gcbo))')];
    
    d=[d,uimenu(c1,'Label','Lost Up','Tag','uifilterauto',...
        'Callback','gui(''uifilter_Callback'',gcbo,14,guidata(gcbo))')];
    d=[d,uimenu(c1,'Label','Lost Dn','Tag','uifilterauto',...
        'Callback','gui(''uifilter_Callback'',gcbo,15,guidata(gcbo))')];

    d=[d,uimenu(c2,'Label','Sort Labels','Tag','uisortplain',...
        'Callback','gui(''uifilter_Callback'',gcbo,16,guidata(gcbo))')];
    
        
    data.uimenu2=d;
        
    
    
    handles.edtlbxeventsold=handles.edtlbxevents;
    handles.fig2=a;
    handles.edtlbxevents=b;        
    
    data.fig2=handles.fig2;
    data.fmain=handles.fmain;
    
    str=get(handles.edtlbxeventsold,'string');
    set(handles.edtlbxevents,'string',str);
    a=get(handles.edtlbxeventsold,'UserData');
    set(handles.edtlbxevents,'UserData',a);    
    a=get(handles.edtlbxeventsold,'Value');
    set(handles.edtlbxevents,'Value',a);
    
    guidata(data.fig2,handles);
    if(~isempty(data.fig1)) guidata(data.fig1,handles); end
    
    % reload content
    reload([],handles);
else
    if(isempty(data.fig2)) return; end
    str=get(handles.edtlbxevents,'string');
    a=get(handles.edtlbxevents,'UserData');    
    b=get(handles.edtlbxevents,'Value');
    
    handles.edtlbxevents=handles.edtlbxeventsold;
    set(handles.edtlbxevents,'string',str);
    set(handles.edtlbxevents,'UserData',a);        
    set(handles.edtlbxevents,'Value',b);
    
    delete(handles.fig2);
    handles.fig2=[];
    data.fig2=[];
    data.uimenu2=[];
    if(~isempty(data.fig1)) guidata(data.fig1,handles); end    
end

guidata(hObject,handles);

% ######################################################################
function uifilter_Callback(hObject, eventdata, handles)
global data

if(isempty(eventdata)) return; end
if(~isfield(data,'uimenu2') || isempty(data.uimenu2)) return; end

entry=data.uimenu2(eventdata);
state=strcmp(get(entry,'Checked'),'on');
if(state) newstate='off'; else newstate='on'; end
set(entry,'Checked',newstate);

% this number should be tracked in the above with
% By Slice filter flag
if((eventdata==8)  & ~state)
    data.uislc=round(get(handles.navsldbrowse,'Value')); 
    set(handles.edtbtnext,'Enable','On');
    set(handles.edtbtprev,'Enable','On');
elseif(eventdata==8)
    set(handles.edtbtnext,'Enable','Off');
    set(handles.edtbtprev,'Enable','Off');
    data.uislc=[];
end

% dependencies
if(~state)
    switch(eventdata)
        case {2,4,14}
            set(data.uimenu2(eventdata+1),'Checked','off');
        case {3,5,15}
            set(data.uimenu2(eventdata-1),'Checked','off');
    end
end

if(eventdata==11 & ~state)
    set(data.uimenu2(16),'Checked','off');
elseif(eventdata==16 & ~state)
    set(data.uimenu2(11),'Checked','off');
end

reload([],handles);

utdgetidx(handles);
utdshow(handles);
        
% ######################################################################
function edtrbtquickdraw_Callback(hObject, eventdata, handles)
% quick draw function
global data

set(handles.edtcklink,'Value',0); data.tmp{4}=[];
set(handles.edtckdelete,'Value',0);
set(handles.edtrbtdraw,'Value',0);
set(handles.edtrbtdelete,'Value',0);
set(handles.edtppcorrect,'Value',1);

% #####################################################################
function edtrbtredraw_Callback(hObject, eventdata, handles)
% saver editing mode - only redraw when this rbutton reclicked
if(get(hObject,'Value'))  utdshow(handles); end

% ######################################################################
function edtrbtthmode_Callback(hObject, eventdata, handles)

% ######################################################################
function edtsavenote(handles,flag)
% function is responsible for saving notes in proof.notes when focus had
% shifted and there is unwritten note left
global proof data debug proof

% get what is the logid in question
noteid=data.edlogid;

% see what is the state of conclusion boxes
istat=get(handles.edtckcorrected,'Value') + ...
             2*get(handles.edtcknotsure,'Value') + ...
              4*get(handles.edtckother,'Value');                        

if(isempty(noteid) || (noteid(1)<=0) || (istat==0)) return; end

% check that the note is saved in proper substack
if(proof.first~=data.edfirst)
    set(handles.navtxtinfo,'String',...
        sprintf('note is for wrong substack, need %i',data.edfirst));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end
    
              
% make another note record and store it
data.ednotes=data.ednotes+1;
proof.notes(data.ednotes).tag=noteid(1);
proof.notes(data.ednotes).slice=noteid(2);
proof.notes(data.ednotes).concl=istat;
proof.notes(data.ednotes).stats=[];
proof.notes(data.ednotes).note=get(handles.edtedtnotes,'String');  

set(handles.navtxtinfo,'String',sprintf('note saved for %i',noteid(1)));
set(handles.navtxtinfo,'fontweight','normal');
pause(0.25);

% don't need to do following if in the wrong substack
if(data.modefirst==proof.first)

    % update tracking
    idd=find(proof.pmap==proof.pmap(noteid(1)+1));
    data.smap(idd)=data.cslc;

    % store timing info
    data.ctime(1)=toc;
    if(data.ctime(2)>0) 
        data.ctime(2)=data.ctime(1)-data.ctime(2); 
    end
    x=proof.ttime(noteid(1),:);
    x([1,2,4:7])=x([1,2,4:7])+data.ctime([1,2,4:7]);
    if(data.ctime(3)>0) x(3)=data.ctime(3); end
    proof.ttime(noteid(1),:)=x; 

    % see if can find in todo list
    % obtain position of the selected fragment in edtlbxevents
    gPos=find(proof.pmap(data.tmp{7}(:,1)+1)==proof.pmap(noteid(1)+1));
    data.tmp{8}(gPos,3)=1;
    
    reload([],handles);
    
    % increase corrections counter
    data.bkcount=data.bkcount+1;    
end

% reset note holders
data.edlogid=[];
data.edconcl=0;
data.edstats=[];
set(handles.edtedtnotes,'String','');

% clean undo if necessary
if(data.pref.undocleanup) data.undo=[]; end

% reset conclusion checkboxes to defaults
set(handles.edtcknotsure,'Value',0);
set(handles.edtckother,'Value',0);
set(handles.edtpptype,'Value',1);

% this to maintain "corrected" selection when cycling through to-do list
if(nargin==1 || flag~=2) set(handles.edtckcorrected,'Value',0); end


% indicate that major processes need to be rebuilt
data.vXupdt=1;


% ######################################################################
function reload(gPos,handles)
% function rolls to-do list to proper position
global data proof
    
if(data.modefirst~=proof.first)
    % deny reloading
    set(handles.navtxtinfo,'String',...
        sprintf('Denied proof list reload in wrong substack (correct %i)',...
         data.modefirst));
    set(handles.navtxtinfo,'fontweight','normal');
    return;
end
    

if(isempty(data.tmp{2}))     
    set(handles.edtlbxevents,'Value',1);
    set(handles.edtlbxevents,'string','nothing to show'); 
    
    % set position vars
    data.ind=[];
    data.lpos=[];
    data.tmp{7}=[];
    
    % sets title-bar
    if(isfield(data,'fig2') && ~isempty(data.fig2))
        set(data.fig2,'name','proofing list');
    end
    
    return;
end

% CHECKED OBJS so far
coms=[proof.notes(1:data.ednotes).tag];

% create internal representation list, if nothing found
if(isempty(data.tmp{7}))
    data.tmp{7}=zeros(length(data.tmp{2}),4);
    data.tmp{7}(:,1)=data.tmp{2};
    % locations of data.major in extended unique list
    [ltmp,idx]=ismember(data.major,data.tmp{2});
    data.tmp{7}(idx(ltmp),2)=data.morder(ltmp);
    data.tmp{7}(idx(ltmp),3)=data.vorder(ltmp);
    data.tmp{7}(idx(ltmp),4)=data.zmajor(ltmp);
    % data.tmp{7} :: obj ID  |  obj ordering  |   obj volume | z-loc ::
    
    % this holds filter/marking options
    data.tmp{8}=false(length(data.tmp{2}),4);
    % data.tmp{8} :: obj filter | obj NS | obj checked | obj shown ::

    % reset initial position
    data.lpos=[];
    data.ind=[];
end

% obtain current position
if(isfield(data,'lpos') && ~isempty(data.lpos))
    % identify the page
    page=max(1,ceil(data.lpos/data.pref.todolength));
    offset=data.lpos-(page-1)*data.pref.todolength;

    nmin=(page-1)*data.pref.todolength+1;
    nmax=min(page*data.pref.todolength,length(data.ind));
    
    % if no in-page position provided, take current
    if(isempty(gPos))
        gPos=get(handles.edtlbxevents,'Value'); 
        % offset for 'prev' entry
        if(nmin>1) gPos=gPos-1; end
    else
        if(nmin>1) gPos=gPos-1; end
    end

    lPos=max(1,min(length(data.ind),nmin+gPos-1));
    glpos=data.ind(lPos);
else
    glpos=[];
end


% obtain filtering status
ufilt=[]; idk=[1:10,14,15];
if(isfield(data,'uimenu2') && ~isempty(data.uimenu2))
    for k=idk ufilt=[ufilt,strcmp(get(data.uimenu2(k),'Checked'),'on')]; end
else ufilt=[1,zeros(1,11)]; end

% obtain sorting status
if(isfield(data,'uimenu2') && ~isempty(data.uimenu2))
    if(strcmp(get(data.uimenu2(11),'Checked'),'on')) usort=1; 
    elseif(strcmp(get(data.uimenu2(16),'Checked'),'on')) usort=2;
    else usort=0; end
    usort=uint8(usort);
else usort=1; end

% obtain skipping status
if(isfield(data,'uimenu2') && ~isempty(data.uimenu2))
    uskip=uint8(strcmp(get(data.uimenu2(12),'Checked'),'on'));
else uskip=get(handles.edtrbtthmode,'Value'); end

% obtain p-map status
if(isfield(data,'uimenu2') && ~isempty(data.uimenu2))
    umap=uint8(strcmp(get(data.uimenu2(13),'Checked'),'on'));
else umap=0; end

% form filtering string
ltmp=true(1,size(data.tmp{8},1));

% apply pmap-transform
if(umap) 
    ltmp=ismember(data.tmp{7}(:,1),proof.pmap(data.tmp{7}(:,1)+1))'; 
end

% major filter
if(ufilt(1)) 
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
    idxtgt=data.major; if(umap) idxtgt=proof.pmap(idxtgt+1); end        
    ltmp=ltmp & ismember(idx,idxtgt); 
end

% all-checked filter
if(ufilt(2) | ufilt(3)) 
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
    idxtgt=coms; if(umap) idxtgt=proof.pmap(idxtgt+1); end        
    if(ufilt(3))
        ltmp=ltmp & ~ismember(idx,idxtgt); 
    elseif(ufilt(2))
        ltmp=ltmp & ismember(idx,idxtgt);
    end
end

% all-typed filter
if(ufilt(4) | ufilt(5)) 
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
    idxtgt=find([proof.tmap]>1)-1; if(umap) idxtgt=proof.pmap(idxtgt+1); end        
    if(ufilt(5))
        ltmp=ltmp & ~ismember(idx,idxtgt); 
    elseif(ufilt(4))
        ltmp=ltmp & ismember(idx,idxtgt);
    end
end

% NS filter

idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
idxtgt=coms(ismember([proof.notes(1:data.ednotes).concl],data.nslist));
idxtgt=idxtgt(:)';
if(umap) idxtgt=proof.pmap(idxtgt+1); end
data.tmp{8}(:,2)=ismember(idx,idxtgt);

if(ufilt(6)) 
    ltmp=ltmp & data.tmp{8}(:,2)';
end

% OTHER filter
if(ufilt(7)) 
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
    idxtgt=coms(ismember([proof.notes(1:data.ednotes).concl],data.otlist));
    if(umap) idxtgt=proof.pmap(idxtgt+1); end
    ltmp=ltmp & ismember(idx,idxtgt); 
end

if(ufilt(11))
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end    
    idxtgt=find(data.objk(:,2)>=data.smax);
    if(umap) idxtgt=proof.pmap(idxtgt+1); end
    ltmp=ltmp & ~ismember(idx,idxtgt); 
end

if(ufilt(12))
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end    
    idxtgt=find(data.objk(:,1)<=1);
    if(umap) idxtgt=proof.pmap(idxtgt+1); end
    ltmp=ltmp & ~ismember(idx,idxtgt); 
end    

% SLICE filter
if(ufilt(8)) 
    % find out where we are
    if(isfield(data,'uislc') && ~isempty(data.uislc))
        kk=data.uislc;
    else
        kk=round(get(handles.navsldbrowse,'Value'));
    end
    
    idx=1:data.ednotes;
    if(ufilt(6))
        idx=idx(ismember([proof.notes(idx).concl],data.nslist));
    end
    if(ufilt(7))
        idx=idx(ismember([proof.notes(idx).concl],data.otlist));
    end
    idx=idx(ismember([proof.notes(idx).slice],kk));    
    
    idxtgt=coms(idx);
    idxtgt=union(idxtgt,data.tmp{7}(data.tmp{7}(:,4)==kk,1));
    if(umap) idxtgt=proof.pmap(idxtgt+1); end
    
    if(ufilt(11))
        idxtgt=find(data.objk(:,2)>data.uislc);
        if(umap) idxtgt=proof.pmap(idxtgt+1); end
    end
    if(ufilt(12))
        idxtgt=find(data.objk(:,1)<data.uislc);
        if(umap) idxtgt=proof.pmap(idxtgt+1); end
    end    
    
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end    
    if(ufilt(11) | ufilt(12))
        ltmp=ltmp & ~ismember(idx,idxtgt);
    else
        ltmp=ltmp & ismember(idx,idxtgt);
    end
    
    if(ufilt(11))
        idxtgt=find(data.objk(:,2)>=data.uislc);
        if(umap) idxtgt=proof.pmap(idxtgt+1); end
    end
    if(ufilt(12))
        idxtgt=find(data.objk(:,1)<=data.uislc);
        if(umap) idxtgt=proof.pmap(idxtgt+1); end
    end
    
    ltmp=ltmp & ismember(idx,idxtgt);
end

% AUTO filter
if(ufilt(10))
    idx=find(ismember([proof.notes(1:data.ednotes).concl],data.nslist));
    ltype=false(1,length(idx));
    for k=1:length(idx)       
        note=proof.notes(idx(k)).note;        
        str={};
        if(~iscell(note))
            for i=1:size(note,1)
                str{end+1}=note(i,:); 
            end
        else
            str=note;
        end
        for i=1:length(str)
            ltype(k)=ltype(k) | ~isempty(strfind(str{i},'AUTO'));
        end
    end
    idxtgt=coms(idx(ltype)); if(umap) idxtgt=proof.pmap(idxtgt+1); end
    
    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
    ltmp=ltmp & ismember(idx,idxtgt); 
end

% TYPE filter
if(ufilt(9))
    % /* time consuming */
    
    % all of this is to increase speed
    ind=data.tmp{7}(ltmp,1);
    ltype=ismember(proof.pmap,proof.pmap(ind+1));
    idx=proof.pmap(ltype);
    idx1=proof.tmap(ltype);
    idx2=idx;
    
    remap=1:max(idx);
    remap(unique(idx))=1:length(unique(idx));
    idx=remap(idx);    
        
    % identify global assignment types
    ss=regionprops(idx,'Area','PixelIdxList');
    for k=find([ss.Area])
        itype=unique(idx1(ss(k).PixelIdxList));
        % primary type overrides where present        
        itype1=itype(ismember(itype,[2,3,4,5]));        
        if(~isempty(itype1)) itype=itype1; end
        ss(k).type=itype;
        ss(k).itype=~isempty(find(ismember(itype,data.pref.objdskt)));
    end
    
    if(~isempty(ss))
        idxtgt=unique(idx2);
        idxtgt=idxtgt([ss.itype]);

        idxtgt=union(idxtgt,data.tmp{7}(...
            ismember(proof.tmap(data.tmp{7}(:,1)+1),data.pref.objdskt),1));
        idxtgt=idxtgt(:)';
        if(umap) idxtgt=proof.pmap(idxtgt+1); end
    else
        idxtgt=[];
    end

    idx=data.tmp{7}(:,1)'; if(umap) idx=proof.pmap(idx+1); end
    ltmp=ltmp & ismember(idx,idxtgt); 
end

% filtered options
data.tmp{8}(:,1)=ltmp;
data.ind=find(ltmp);
lcoms=data.tmp{8}(ltmp,2);
lchk=data.tmp{8}(ltmp,3);

% sort if requested
switch(usort)
    case 1
        idx=data.tmp{7}(ltmp,2);
        [idx1,idx]=sort(idx,'descend');
        data.ind=data.ind(idx);
        lcoms=lcoms(idx);
        lchk=lchk(idx);
    case 2
        idx=data.tmp{7}(ltmp,1);
        [idx1,idx]=sort(idx,'descend');
        data.ind=data.ind(idx);
        lcoms=lcoms(idx);
        lchk=lchk(idx);        
end

ind=data.tmp{7}(data.ind,1);
if(isempty(ind))     
    set(handles.edtlbxevents,'Value',1);    
    set(handles.edtlbxevents,'string','nothing to show');
    
    % set position vars
    data.ind=[];
    data.lpos=[];

    % sets title-bar
    if(isfield(data,'fig2') && ~isempty(data.fig2))
        vorder=data.tmp{7}(:,3);

        % estimate processed volume
        v=sum(vorder(ismember(proof.pmap(data.tmp{7}(:,1)+1),coms)));
        v=round(100*v/sum(vorder));

        % largest order. val. below current pos
        s=sprintf('*%i* N:%i M:0 R:%i%%',0,length(unique(coms)),0,v);
        set(data.fig2,'name',s);
    end
    
    return;
end

% identify current position
if(isempty(glpos)) glpos=data.ind(1); end

% do this b/s obj may not be in the list anymore
if(~ismember(glpos,data.ind))
    if(usort)
        idx=data.tmp{7}(ltmp,2);        
        idx=abs(idx-data.tmp{7}(glpos,2)); 
        [a,mm]=min(idx);
    else
        idx=abs(ind-data.tmp{7}(glpos,1)); [a,mm]=min(idx1);
    end
    glpos=data.ind(mm);
end

% find out on which page this is located & local position in filtered list
% local position could have changed per filtering
lpos=find(data.ind==glpos);
data.lpos=lpos;

% identify local page
page=max(1,ceil(lpos/data.pref.todolength));
offset=lpos-(page-1)*data.pref.todolength;

nmin=(page-1)*data.pref.todolength+1;
nmax=min(page*data.pref.todolength,length(ind));

% fill out the page
coms=proof.pmap(coms+1);
str=cell(nmax-nmin+1,1);
if(nmin>1) str{1}='   prev    '; di=1; else di=0; end
for i=1:nmax-nmin+1
    nclst=ind(i+nmin-1);
    if(lcoms(i+nmin-1)) s='*comment*'; else s=' '; end
    if(ismember(proof.pmap(nclst+1),coms)) s1='C';  else s1=' '; end
    itype=proof.tmap(nclst+1); itype1=[];
    if(itype==1)
        itype=unique(proof.tmap(proof.pmap==proof.pmap(nclst+1)));
        itype1=unique(itype(ismember(itype,[2,3,4,5])));
        if(~isempty(itype1)) itype=itype1; end
        itype=itype(itype>1);
    end
    if(length(itype1)>1) 
        stype='conflict'; 
    elseif(length(itype)==1)
        stype=get(handles.edtpptype,'string');
        stype=stype{itype};
    else stype=''; end
    
    str{i+di}= sprintf('%s #%i %8s %6i [%i] %s',...
        s1,i+nmin-1,stype,nclst,lchk(i+nmin-1),s);
end
if(nmax<length(ind)) str{end+1}='   next   '; end
set(handles.edtlbxevents,'String',str);
if(nmin>1) offset=offset+1; end
set(handles.edtlbxevents,'Value',offset);

% sets title-bar
if(isfield(data,'fig2') && ~isempty(data.fig2))
    vorder=data.tmp{7}(:,3);

    % estimate processed volume
    v=sum(vorder(ismember(proof.pmap(data.tmp{7}(:,1)+1),coms)));
    v=v/sum(vorder);

    % largest order. val. below current pos
    morder=data.tmp{7}(:,2);
    [ind,idx]=sort(morder,'descend');
    m=ind(find(idx==glpos)+1);
    if(m>1E6)
        s=sprintf('%iM',ceil(m/1e6));
    elseif(m>1E3)
        s=sprintf('%iK',ceil(m/1e3));
    else
        s=sprintf('%i',ceil(m));
    end
    s=sprintf('*%i* N:%i M:%s R:%i%%',length(find(ltmp)),...
        length(unique(coms)),s,round(100*v));
    set(data.fig2,'name',s);
end


% ######################################################################
function edtselection(handles)
% function resets editing controls at new selection of target segment
global proof data debug notes

if(data.curfirst~=proof.first)
    set(handles.navtxtinfo,'string',...
        sprintf('you are in wrong substack, need %i',data.curfirst));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

% see if we are in guided-inspection mode
imode=get(handles.navppselect,'Value');
imode=data.selmodes(imode);
iflag=(imode==6);

% obtain position of the selected fragment in edtlbxevents
gPos=find(proof.pmap(data.tmp{7}(:,1)+1)==proof.pmap(data.nclst+1));

% take special care of what is selected if supervising;
%  more than one record may be present
pmode=get(handles.edtppmode,'Value');
pmode=data.edtmodes(pmode);

% in case there are more than single record
if(length(gPos)>1)
    if(ismember(data.ind(data.lpos),gPos)) 
        gPos=data.ind(data.lpos);
    else gPos=gPos(1); end
end


if(gPos>0) % if there is such...
    data.tmp{18}=data.tmp{7}(gPos,1);
    % see if we need to save anything before going on to another fragment
    edtsavenote(handles,data.pref.dftcor);

    % reset editing controls to their defaults
    set(handles.navppselect,'Value',find(data.selmodes==3));
    set(handles.edtcknotsure,'Value',0);
    set(handles.edtckother,'Value',0);
    set(handles.edtpptype,'Value',1);
    set(handles.edtppcorrect,'Value',1);
%     set(handles.edtlbxevents,'Value',iPos);
    set(handles.edtcklink,'Value',0);   
    set(handles.edtckdelete,'Value',0);       
    set(handles.edtrbtdraw,'Value',0);
    set(handles.edtrbtdelete,'Value',0);
    data.tmp{4}=[];
    % reset bwlabel for quick-draw
    data.tmp{17}=[];    
        
    % use preference for corrected unless in supervising mode
	switch(data.pref.dftcor)
        case 0
            set(handles.edtckcorrected,'Value',0);
        case 1
            set(handles.edtckcorrected,'Value',pmode~=9);
    end
    % use preference for pen size
    set(handles.edtppsize,'Value',data.pref.dftpen);
    % see if quick-link should be by default enabled
    if(data.pref.dftcklink & pmode~=9)
        set(handles.edtcklink,'Value',1);
        set(handles.edtppcorrect,'Value',2);
        data.tmp{4}=data.tmp{7}(gPos,1);
    end
        
    % if we are at the end of the list, disable cycle button
    if(data.lpos==length(data.ind))
        set(handles.edtbtcycle,'enable','off');
    else
        set(handles.edtbtcycle,'enable','on');
    end

    % extract labels associated with the fragments
    % what is the fragment's slice
    k=data.tmp{5};
    % what is the fragment's primary tag
    ntag=0;
    % what is the fragment's outgoing tag
    nclst0=data.tmp{7}(gPos,1);
    % what is the fragments global id??? may be same with outgoing
    nclst=data.tmp{7}(gPos,1);            
    
    % output debug info
    set(handles.dbgedtstats,'String',{''});
    data.logid=dbginfo(ntag,nclst0,nclst,k,handles);
    
    % add comment entries if present
    str=get(handles.dbgedtstats,'String');
    if(~iscell(str)) str={str}; end
    coms=[proof.notes(1:data.ednotes).tag];
    icoms=find(ismember([proof.notes(1:data.ednotes).concl],data.nslist));
    coms=coms(icoms);    
    if(ismember(data.nclst,coms))
        icur=find(coms==data.nclst);
    elseif(ismember(proof.pmap(data.nclst+1),proof.pmap(coms+1)))
        icur=find(proof.pmap(coms+1)==proof.pmap(data.nclst+1));
    else
        icur=[];
    end
    for i=1:length(icur)
        j=icoms(icur(i));
        str{end+1}=sprintf('*** NS PROOF-NOTE #%i CONCLUSION %i SLICE %i***',...
            j,proof.notes(j).concl,proof.notes(j).slice);
        if(isempty(str{1})) str={str{2}}; end
        s=proof.notes(j).note;
        if(isempty(s)) s='none'; end
        str=cat(1,str,s);
    end
    set(handles.dbgedtstats,'String',str);    
%     if(pmode~=9) set(handles.dbgedtstats,'String',str); end
    
%     % if in supervisor mode, output one notes record
%     % SUPERVISOR MODE IS DISABLED
%     if(pmode==9 && ~isempty(data.tmp{11}))
%         irec=data.tmp{11}(gPos);
%         
%         inote=proof.notes(irec);
%         if(ismember(inote.concl,data.clist))
%             set(handles.edtckcorrected,'Value',1);
%         end
%         
%         if(ismember(inote.concl,data.nslist))
%             set(handles.edtcknotsure,'Value',1);
%         end
%         
%         if(ismember(inote.concl,data.otlist))
%             set(handles.edtckother,'Value',1);
%         end
%         
%         if(~isempty(inote.note))
%             if(iscell(inote.note)) str1=inote.note; else str1={inote.note}; end
%         else
%             str1={''};
%         end
%         str1=cat(1,{sprintf('*** NS PROOF-NOTE #%i CONCLUSION %i SLICE %i***',...
%             irec,proof.notes(irec).concl)},proof.notes(irec).slice,str1);
%         
%         set(handles.dbgedtstats,'String',str1);
%         
%         data.edlogid=[];
%     end      
    
    % reset sensitive variables
    data.edstats=[];        % edstats, edlogid and edconcl for notes
    if(get(handles.edtckcorrected,'Value') & (pmode~=9))
        data.edlogid=[nclst,k];
        data.edconcl=1;
    else
        data.edlogid=[];
        data.edconcl=0;
    end
        
    % set tracking to this fragment if user will wish to track it
    % as a cluster independently
    data.ids=proof.pmap(1+nclst);
    data.nclst=nclst;
    set(handles.navedtselection,'string',num2str(nclst));

    id=nclst;
    itype=proof.tmap(id+1);
    if(itype==1)
        itype=proof.tmap(proof.pmap==proof.pmap(id+1));
        itype=intersect(itype,[2,3,4,5]);
        if(length(itype)>1) itype=1; 
        elseif(length(itype)==0) itype=proof.tmap(id+1);  end
    end
    set(handles.edtpptype,'value',itype);
    
    s=get(handles.edtpptype,'String');
    s=sprintf('selected %s%i',s{proof.tmap(nclst+1)}(1),nclst);
    set(handles.navtxtinfo,'string',s);
    set(handles.navtxtinfo,'fontweight','normal');

    % tell user when this process was last time saved
    str=get(handles.navtxtinfo,'string');
    if(iscell(str)) str=str{1}; end
    nn=ismember(proof.pmap(1+nclst),...
        proof.pmap(1+[proof.notes(1:data.ednotes).tag]));
    nn=max(double(nn),data.smap(1+nclst));
    if(nn>0)
        set(handles.navtxtinfo,'String',[str,': last edt ',num2str(nn)]);
        set(handles.navtxtinfo,'fontweight','bold');
    end
    
    % movein on the fragment, only if we actually in guided mode
    if(iflag)
        data.v2zoom=[];
        utdshow(handles);
        utdzoomin(handles);
    end   
        
    % reset timing var
    data.ctime=zeros(1,7);
    % start timer
    tic
    % update editing starting time
    if(get(handles.edtckcorrected,'value')) data.ctime(2)=1e-6; end
end

% ######################################################################
function edtundo(handles)
% this function handles undo capacity
global cat proof data

i=length(data.undo);
if(i==0)
    set(handles.navtxtinfo,'String','undo stack is empty');
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end
undo=data.undo(i);

% check we are in correct substack 
if(undo.proof~=proof.first)
    set(handles.navtxtinfo,'String',sprintf('you are in wrong substack, need %i',undo.proof));
    set(handles.navtxtinfo,'fontweight','bold');
    return;
end

switch(undo.action)
    case {1,2,3,4}
        tmp=cat{undo.slc};
        tmp(undo.data{2})=undo.data{1};
        cat{undo.slc}=tmp;
        data.undo=data.undo(1:end-1);
    case {5,6}
        proof.pmap(undo.data{2})=undo.data{1};
        data.undo=data.undo(1:end-1);
        % tell GUI that significant object may need to be recalculated
        data.vXupdt=1;        
end
set(handles.navtxtinfo,'String','undo');
set(handles.navtxtinfo,'fontweight','normal');
        
utdgetidx(handles);
utdgetimg(handles);

% redraw to show corrections made
if(get(handles.edtrbtredraw,'Value')) utdshow(handles); end



% #######################################################################
% ############          AUXILIARY

% ######################################################################
function mainKeyDownFcn(hObject,handles)
% key-shortcuts handler
global data

% what was clicked
s=double(get(handles.fmain,'CurrentCharacter'));
if(isempty(s)) return; end

% were we are right now
cntr=ceil([data.zoom(1)+data.zoom(2),data.zoom(3)+data.zoom(4)]/2);;
DD=[data.zoom(2)-data.zoom(1),data.zoom(4)-data.zoom(3)]+1;
DD=max(5,floor(DD/10));
k=round(get(handles.navsldbrowse,'Value'));

% load up preferences
A=data.pref.keys;
% this is hot-switch to primary object type
if(~ismember(s,A) && (s>='1' & s<='3'))
    set(handles.edtpptype,'value',uint8(s)-47);
    edtpptype_Callback(handles.edtpptype,[],handles);
end

% nongeneric work-classification module
if(~ismember(s,A) && (s>='4' & s<='7'))
    msg={'trace/solve','alignment','vesicles','mitochondria'};
    x=double(s)-51;
    data.ctime(3)=x;
    set(handles.navtxtinfo,'String',['work class ',msg{x}]);
    set(handles.navtxtinfo,'fontweight','normal');
end 
            
switch(s)
    case A(1) % up key
        cntr(1)=max(1,cntr(1)-DD(1));
        utdzoom(cntr,0.5,handles);
    case A(2) % down key
        cntr(1)=min(data.shape(1),cntr(1)+DD(1));
        utdzoom(cntr,0.5,handles);
    case A(3) % left key
        cntr(2)=max(1,cntr(2)-DD(2));
        utdzoom(cntr,0.5,handles);
    case A(4) % right key
        cntr(2)=min(data.shape(2),cntr(2)+DD(2));
        utdzoom(cntr,0.5,handles);
    case A(5) % '-'-key -- zoom in
        utdzoom(cntr,0.25,handles);
    case A(6) % '='-key -- zoom out
        utdzoom(cntr,1.0,handles);
    case A(7) % ']'-key -- move up
        utdmove(k+1,handles);
    case A(8) % '['-key -- move down
        utdmove(k-1,handles);
    case A(9) % 'f'-key -- flip
        navrbtflip_Callback([],[],handles);
    case A(10) % 'c'-key -- cycle
        edtbtcycle_Callback(handles.edtbtcycle,[],handles);        
    case A(11)% 'o'-key -- overlay
        set(handles.navckoverlay,'Value',~get(handles.navckoverlay,'Value'));
        navckoverlay_Callback(handles.navckoverlay,[],handles);
    case A(12)% 'n'-key -- next
        edtbtnext_Callback(handles.edtbtnext,[],handles);
    case A(13)% 'p'-key -- prev
        edtbtprev_Callback(handles.edtbtprev,[],handles);
    case A(14)% 'd'-key -- quick draw
        set(handles.edtrbtquickdraw,'Value',1);
        edtrbtquickdraw_Callback(handles.edtrbtdraw,[],handles);
    case A(15)% 'q'-key -- quick link
        set(handles.edtcklink,'Value',1);
        edtcklink_Callback(handles.edtcklink,[],handles);
    case A(16)% 'r'-key -- redraw
        imode=logical(get(handles.edtrbtredraw,'Value'));
        set(handles.edtrbtredraw,'Value',~imode);
        edtrbtredraw_Callback(handles.edtrbtredraw,[],handles);
    case A(17)% 'z'-key -- undo
        edtundo(handles);
    case A(18)% 's'-key -- quick del
        set(handles.edtckdelete,'value',1);
        edtckmembrane_Callback(handles.edtckdelete,[],handles);
    case A(19)% 'e'-key -- del
        set(handles.edtrbtdelete,'value',1);
        edtrbtdelete_Callback(handles.edtrbtdelete,[],handles);
    case A(20)% 'r'-key -- draw
        set(handles.edtrbtdraw,'value',1);
        edtrbtdraw_Callback(handles.edtrbtdelete,[],handles);
    case A(21)% 'x'-key -- recenter
        utdzoomin(handles);
end

% display(s);


function mkwatershed(handles)
% implements watershed functionality
% this script redraws proof-read segmentation using watershed 
%  from cleaned fragments;
% Input:
%   al & mskE & proof & debug or major"=dbggetmajor(debug,proof)"
% executing inside
%     cat"=proof.pmap(imadd(cat,1))"
% Ouput: 
%   wcat stack of finalized fragments
%   ucat cleaned stack of only major fragments
global al cat debug proof data

% constants: [close holes, step back to define seeds, size of holes to close]
H=[25 15 200]*data.pref.edrange/2;
% flag to fix possibly broken 4-borders
fix4=0;
% current slice
kk=round(get(handles.navsldbrowse,'Value'));
ikk=kk-data.cpos+1;

str{1}='executing watershed...';
str{2}='extracting major processes...';
set(handles.dbgedtstats,'String',str);
pause(0.1);

% get major objects
if(data.vXupdt)
    proof.major=dbggetmajor;
    data.vXupdt=0;
end
major=proof.pmap(proof.major+1);

str{3}=sprintf('major processes: %i',length(major));
set(handles.dbgedtstats,'String',str);

% filter to define LOG edges
h=fspecial('gauss',13,2);

% STD search
h1=ones(3,3)/9;

% diameter to fill in holes to define actual area occupied by 
% found fragments (this may be smaller than mskE b/s of fragments
% on the boundary erased in proof-reading)
se1=strel('square',H(1));

% diameter to step back from thus above dilated region to define
% seed for boundary-connected-background
se2=strel('square',H(2));

if(isfield(debug,'ind1') && ~isempty(debug.ind1))
    str=get(handles.dbgedtstats,'String');
    if(~iscell(str)) str={str}; end
    str{end+1}='overwriting debug.ind1';
    set(handles.dbgedtstats,'String',str);
    pause(0.1);
end    

switch(data.pref.watershed)
    case 1 % redraw single section
        idk=kk;
        debug.ind1{kk}=[];
    case 0 % redraw all sections
        idk=1:length(cat);
        debug.ind1=cell(size(cat));        
end
        
for k=idk
    ik=k-data.cpos+1;
    ifget('alcat',ik);    
    
    str=get(handles.dbgedtstats,'String');
    if(~iscell(str)) str={str}; end
    if(strcmp(str{1},'') | strcmp(str{1},'none')) str={}; end
    l=length(str);
    str{l+1}=sprintf('redrawing #%i...',k);
    set(handles.dbgedtstats,'String',str);
    pause(0.1);
    
    % fix possibly broken 1pxl borders
    if(fix4)
        str{end+1}='fixing broken 1pxl borders';
        set(handles.dbgedtstats,'String',str);
        pause(0.1);
        wtmp=cat{ik};
        ltmp=imdilate(wtmp,ones(3,3))~=wtmp;
        wtmp(wtmp==0)=data.mmax+1;
        ltmp=ltmp | (imerode(wtmp,ones(3,3))~=wtmp);
        
        wtmp(ltmp | cat{ik}==0)=0;
    else
        wtmp=cat{ik};
    end

    % obtain mask of only significant objects
    ladd=ismember(proof.pmap(wtmp+1),major);
    
    % define watershed markers and watershed borders
    dtmp=imcomplement(im2double(al{ik})); 
    dtmp=imfilter(dtmp,h);
    
    dtmp=dtmp + (wtmp==0) + (ladd==0);

    ltmp=ladd>0;
    ltmp=imreconstruct(ltmp,wtmp>0);

    ltmp1=imdilate(ltmp,se1);
    wtmp1=bwlabel(~ltmp1);
    stats1=regionprops(wtmp1,'MajorAxisLength');
    ltmp2=ltmp1 | ismember(wtmp1,find([stats1.MajorAxisLength]<H(3)));    
%     ltmp2=imfill(ltmp1,'holes');
    
    ltmp1=~(imerode(ltmp2,se2) | (ltmp2 & ~ltmp1));
    ltmp1=ltmp1 | imreconstruct(ltmp1 & wtmp>0,ltmp1 | wtmp>0);
    
    dtmp=imimposemin(dtmp,ltmp1 | ltmp);
    dtmp=watershed(dtmp,8);

    % renumber watershed regions with their cat-labels    
    stats1=regionprops(dtmp,'PixelIdxList');
    pmap=linspace(0,length(stats1),length(stats1)+1);
    wtmp1=proof.pmap(wtmp+1);
    wtmp1(~ismember(wtmp1,major))=0;
    for l=1:length(stats1)
        tmp1=double(sort(wtmp1(stats1(l).PixelIdxList(:))));
        tmp2=diff([tmp1;max(tmp1)+1]);
        ssum1=diff(find([1;tmp2]));
        sidx=tmp1(tmp2>0);
        ssum=ssum1(sidx>0);
        ccidx=uint32(sidx(sidx>0));

        [grb,pos]=max(ssum);
        if(length(pos)==0) pos=0; else pos=ccidx(pos); end
        pmap(l+1)=pos;
    end    
    dtmp=pmap(dtmp+1);
    
    % merge neighbour watershed fragments together
    se=ones(3,3);
    dtmp=imdilate(dtmp,se);
    ltmp=(dtmp~=imdilate(dtmp,se)) | (dtmp~=imerode(dtmp,se));
    dtmp(ltmp)=0;
    if(isfield(debug,'mskE') && ~isempty(debug.mskE) && ~isempty(debug.mskE{k}))
        dtmp(debug.mskE{k})=0;
    end
    
    if(data.mmax<2^16) debug.ind1{k}=uint16(dtmp); else debug.ind1{k}=dtmp; end
    
    % save-to-file
    if(isfield(debug,'prefix') && ~isempty(debug.prefix))
        sname=sprintf('%s.%.3i.mat',debug.prefix,k);
        wcat=debug.ind1{k};
        save(sname,'-append','wcat');
        debug.ind1{k}=[];
        wcat=[];
    end
end


% set-up display mode
str=get(handles.dbgedtstats,'String');
if(~iscell(str)) str={str}; end
if(strcmp(str{1},'') | strcmp(str{1},'none')) str={}; end
l=length(str);
str{l+1}='finished';
set(handles.dbgedtstats,'String',str);
pause(0.1);

set(handles.navckoverlay,'Value',1);
data.tmp{1}=6;
set(handles.navppimgmode,'Value',find(data.dismodes==16));
set(handles.navrbtwtshdcntrs,'Value',1);

data.v1zoom=[];
data.v2zoom=[];
utdshow(handles);












% ######################################################################
function savebackup(mode,handles)
% function designed to load/save backup
global proof cat data debug

global flMonitor

if(~isempty(flMonitor) && flMonitor && mode ==0)
%     set(handles.navtxtinfo,'String',...
%                     'save disabled when monitoring, set fmonitor=0');
%     set(handles.navtxtinfo,'fontweight','bold');
    return;
end

% check that the note is saved in proper substack
if(proof.first~=data.edfirst)
    set(handles.navtxtinfo,'String',...
        sprintf('can''t do - wrong substack, need %i',data.edfirst));
    set(handles.navtxtinfo,'fontweight','bold');
end

% record state of gui
data.internal.dismode=get(handles.navppimgmode,'Value');
data.internal.edtmode=get(handles.edtppmode,'Value');
data.internal.selmode=get(handles.navppselect,'Value');

data.internal.thmode=get(handles.edtrbtthmode,'Value');

data.internal.edtlist=get(handles.edtlbxevents,'String');
data.internal.edtpos=get(handles.edtlbxevents,'Value');
data.internal.edtnote=get(handles.edtedtnotes,'String');

data.internal.imgpos=[get(handles.navrbtmap,'Value'),...
    get(handles.navckoverlay,'Value'),get(handles.navckoverlap,'Value')];
data.internal.selection=get(handles.navedtselection,'String');
data.internal.seltype=get(handles.navrbtselstyle,'Value');

data.internal.slice=get(handles.navsldbrowse,'Value');
data.internal.mix=get(handles.navsldmixer,'Value');

data.internal.nextstate=get(handles.edtbtnext,'enable');
data.internal.prevstate=get(handles.edtbtprev,'enable');
data.internal.cyclestate=get(handles.edtbtcycle,'enable');

data.internal.todopage=get(handles.edtlbxevents,'userdata');

% save full backup every 10 saves
if(mode==0) data.bksaves=data.bksaves+1; end
if(mode==2) data.bksaves=0; end

if(data.bksaves>10 & mode==0) mode=2; data.bksaves=0; end

% current section
data.sectn=get(handles.navsldbrowse,'Value');
switch(mode)
    case 0
        data.saved=1;
        data.bkcount=0;        
        set(handles.navtxtinfo,'string','saving backup pls wait...');
        set(handles.navtxtinfo,'fontweight','bold');
        pause(0.1);
        
        s=sprintf('gui.backup%.3i.mat',proof.first);
        if(exist(s))
            sX=sprintf('gui.backup%.3iX.mat',proof.first);
            copyfile(s,sX);
        end
        % do this to remove drawing saves from short-backup
        tmp=data.tmp{20}; data.tmp{20}={};
        save(s,'proof','data'); data.tmp{20}=tmp;        
        data.undo=[];
        
        set(handles.navtxtinfo,'string','saved'); 
        set(handles.navtxtinfo,'fontweight','normal');
    case 1
        set(handles.navtxtinfo,'string','loading backup pls wait...');
        set(handles.navtxtinfo,'fontweight','bold');
        pause(0.25);
        
        
        s=sprintf('gui.backup%.3i.mat',proof.first);
        load(s);
        
        wcat=[]; ind=[];        
        s=sprintf('gui.data%.3i.mat',proof.first);
        load(s,'ind','wcat');
        data.tmp{20}=cell(size(cat));
        for k=1:length(ind)
            if(~isempty(ind{k}))
                data.tmp{20}{k}=ind{k};
                if(~isempty(cat{k})) cat{k}(ind{k}{1})=ind{k}{2}; end
            end                
        end
        if(~isempty(wcat)) debug.ind1=wcat; end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPATIBILITY SWITCHES
        if(~isfield(data.pref,'keydscrp'))
            data.pref.keydscrp={'up';'down';'left';'right';'zoom in';...
                'zoom out';'move up';'move down';'flip';'cycle';...
                'overlay';'next';'previous';'quick draw';'quick link';...
                'redraw';'undo';'quick del';'delete';'draw';'recenter'};
            data.pref.keys=[data.pref.keys,120];
        end
        if(~isfield(data,'seppos'))
            data.seppos=get(handles.fmain,'Position'); 
        end
        if(~isfield(data,'ctime')) data.ctime=zeros(1,7); end
        if(~isfield(proof,'ttime')) 
            proof.ttime=zeros(length(proof.pmap),7); 
        end
        if(~isfield(data,'vorder')) data.vorder=data.morder; end
        if(~isfield(data,'zmajor')) 
            data.zmajor=zeros(size(data.major),'uint16'); 
        end

        % create internal representation list, if nothing found
        if(isempty(data.tmp{7}) && ~isempty(data.tmp{2}))
            fprintf('WARNING: v7.0 or prior backup is being loaded;\n');
            fprintf('         position in the proofing list will be reset!\n');
            data.tmp{7}=zeros(length(data.tmp{2}),4);
            data.tmp{7}(:,1)=data.tmp{2};
            % locations of data.major in extended unique list
            [ltmp,idx]=ismember(data.major,data.tmp{2});
            data.tmp{7}(idx(ltmp),2)=data.morder(ltmp);
            data.tmp{7}(idx(ltmp),3)=data.vorder(ltmp);
            data.tmp{7}(idx(ltmp),4)=data.zmajor(ltmp);
            % data.tmp{7} :: obj ID  |  obj ordering  |   obj volume | z-loc ::

            % this holds filter/marking options
            data.tmp{8}=false(length(data.tmp{2}),4);
            % data.tmp{8} :: obj filter | obj NS | obj checked | obj shown ::
            if(proof.v==7) data.tmp{8}(:,3)=(data.tmp{9}>0); end

            % reset initial position
            data.lpos=[];
            data.ind=[];
        end
        
        if(~isfield(data,'uislc')) data.uislc=[]; end
        
        if(~isfield(data,'shadow')) data.shadow=0; end
        
        if(~isfield(data,'altcol')) data.altcol=0; end
        
        if(~isfield(data,'modefirst')) data.modefirst=proof.first; end
        
        % update array length array
        data.mmax=max(data.mmax,length(proof.pmap));
        data.mmax=max(data.mmax,double(max(proof.pmap)));        
        data.cmap=[];
        
        % compatibility switch for shading convention
        if(data.pref.shade(1)<1) data.pref.shade(1)=1/data.pref.shade(1); end
        if(data.pref.shade(2)<1) data.pref.shade(2)=1.1; end
        % COMPATIBILITY SWITCHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % restore state
        if(~isempty(data.internal))
            set(handles.navppimgmode,'Value',data.internal.dismode);
            set(handles.edtppmode,'Value',data.internal.edtmode);
            set(handles.navppselect,'Value',data.internal.selmode);

            set(handles.edtlbxevents,'String',data.internal.edtlist);
            set(handles.edtlbxevents,'Value',data.internal.edtpos);
            set(handles.edtedtnotes,'String',data.internal.edtnote);

            set(handles.navrbtmap,'Value',data.internal.imgpos(1));
            set(handles.navckoverlay,'Value',data.internal.imgpos(2));
            set(handles.navckoverlap,'Value',data.internal.imgpos(3));
            set(handles.navedtselection,'String',data.internal.selection);
            set(handles.navrbtselstyle,'Value',data.internal.seltype);
            
            set(handles.edtbtnext,'enable',data.internal.nextstate)
            set(handles.edtbtprev,'enable',data.internal.prevstate);
            set(handles.edtbtcycle,'enable',data.internal.cyclestate);

            set(handles.edtrbtthmode,'Value',data.internal.thmode);
            set(handles.edtlbxevents,'userdata',data.internal.todopage);

            set(handles.navsldmixer,'Value',data.internal.mix);
            utdmove(data.internal.slice,handles);            
        end        
        data.tmp{4}=[];
        
        % reset editing controls
        set(handles.edtckcorrected,'Value',0);
        set(handles.edtcknotsure,'Value',0);
        set(handles.edtckother,'Value',0);
        
        % reset references to outdated handles
        set(handles.edtrbtproofing,'Value',0);
        if(~isempty(data.fig2) &...
                (~isfield(handles,'fig2') || isempty(handles.fig2)))
            data.fig2=[];
        end
        edtrbtproofing_Callback(handles.edtrbtproofing,[],handles);
        
        set(handles.navrbtconsole,'Value',0);
        if(~isempty(data.fig1) &...
                (~isfield(handles,'fig1') || isempty(handles.fig1)))
            data.fig1=[];
        end        
        navrbtconsole_Callback(handles.navrbtconsole,[],handles);
        
        data.uimenu2=[];
        
        data.cntxt=[];
        
        set(handles.edtcklink,'Value',0);
        set(handles.edtckdelete,'Value',0);
        set(handles.edtrbtquickdraw,'Value',0);
        set(handles.edtrbtdraw,'Value',0);
        set(handles.edtrbtdelete,'Value',0);
        
        set(handles.edtppcorrect,'Value',1);
        
        data.bkcount=0;        
        set(handles.navtxtinfo,'string','loaded');   
        set(handles.navtxtinfo,'fontweight','normal');
        
        data.shadow=0;
        
        % reset timer counter
        tic;
    case 2
        data.saved=1;
        data.bkcount=0;        
        set(handles.navtxtinfo,'string','saving backup pls wait...');
        set(handles.navtxtinfo,'fontweight','bold');
        pause(0.1);
        s=sprintf('gui.backup%.3i.mat',proof.first);
        if(exist(s))
            sx=sprintf('gui.backup%.3iX.mat',proof.first);
            copyfile(s,sx);
        end
        % save state
        tmp=data.tmp{20}; data.tmp{20}={};
        save(s,'proof','data');
        data.tmp{20}=tmp;
        
        s=sprintf('gui.data%.3i.mat',proof.first);
        if(exist(s))
            sx=sprintf('gui.data%.3iX.mat',proof.first);
            copyfile(s,sx);
        end                

        % save current proofing status
        sbuf=sprintf('gui.buff%.3i.mat',proof.first);

        ind={};
        if(isempty(data.tmp{20})) data.tmp{20}=cell(size(cat));   end
        for k=1:length(cat) ind{k}=data.tmp{20}{k};   end
        if(isfield(debug,'ind1') && ~isempty(debug.ind1))
            wcat=debug.ind1;
        else
            wcat={};
        end
        
        % disallow saving debug.ind1 if communicated via HD
        if(isfield(debug,'prefix') && ~isempty(debug.prefix)) wcat={}; end
        
        save(s,'ind','proof','wcat');
        save(sbuf,'-v6','ind','proof','wcat');

        data.undo=[];
        set(handles.navtxtinfo,'string','saved'); 
        set(handles.navtxtinfo,'fontweight','normal');        
end



% ======================================================================
%                   MATLAB CONTROLS CREATION
% ======================================================================

function edtppmode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function navsldbrowse_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function navsldmixer_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function navedtselection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dbgedtstats_CreateFcn(hObject, eventdata, handles)
 if ispc && isequal(get(hObject,'BackgroundColor'), ...
         get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
 end

function navppimgmode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function navppselect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dbgpplog_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dbgppevents_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtppcorrect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edtedtnotes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtppevents_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtppsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtpptype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function navedtgoto_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtlbxevents_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
         get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% ###################### ADDONS #########################
% --- Executes on button press in navrbtcross.
function navrbtcross_Callback(hObject, eventdata, handles)

% #######################################################################
function uimpref_Callback(hObject, eventdata, handles)
% calls preference editor
global data

data.fmain=handles.fmain;
if(~isfield(data,'cntxt') || isempty(data.cntxt)) 
    data.cntxt=get(handles.axmain,'UIContextMenu'); 
end

pref;

% #######################################################################
function uimclearind1_Callback(hObject, eventdata, handles)
global data debug
k=round(get(handles.navsldbrowse,'Value')); % current section id
if(data.pref.watershed) debug.ind1{k}=[]; else debug.ind1={}; end
utdgetimg(handles);
utdshow(handles);

% #######################################################################
function navrbtwtshdcntrs_Callback(hObject, eventdata, handles)
global data

if(get(hObject,'value'))
    % set-up display mode
    set(handles.navckoverlay,'Value',1);
    data.tmp{1}=6;
    set(handles.navppimgmode,'Value',find(data.dismodes==16));

    data.v1zoom=[];
    data.v2zoom=[];
    utdshow(handles);
else
    % set-up display mode
    set(handles.navckoverlay,'Value',0);
    data.tmp{1}=[];
    set(handles.navppimgmode,'Value',find(data.dismodes==1));

    data.v1zoom=[];
    data.v2zoom=[];
    utdshow(handles);
    set(hObject,'value',0);    
end


% --- Executes when fmain is resized.
function fmain_ResizeFcn(hObject, eventdata, handles)
% this is manual resize function
global data

if(~isfield(data,'seppos')) return; end

cclist=get(handles.fmain,'Children');

% width of the main window
w=get(handles.fmain,'Position');
DX=w(3)-data.seppos(3);
DY=w(4)-data.seppos(4);

% axmain needs separate handler
if(~isfield(data,'axmainsize') || isempty(data.axmainsize))
    data.axmainsize=get(handles.axmain,'Position');
end

% controls are tied to left edge
for i=1:length(cclist)
    % only affect "non-proportional" controls
    if(isprop(cclist(i),'Units')) u=get(cclist(i),'Units'); else u=''; end
    if(strcmp(u,'pixels'))
        % get original control specs
        ud=get(cclist(i),'UserData');
        up=get(cclist(i),'Position');
        if(isempty(ud))
            up=get(cclist(i),'Position');
            set(cclist(i),'UserData',up);
            ud=up;
        end
        up(1)=ud(1)+DX;
        up(2)=ud(2)+DY;
        
        set(cclist(i),'Position',up);
    end
end

% main axis scales with the window to fill in the area
ud=data.axmainsize;
up=ud; 
up(3)=up(3)+DX; up(4)=up(4)+DY;
set(handles.axmain,'Position',up);

% adjust slider bars
ud1=get(handles.navsldbrowse,'UserData');
ud2=get(handles.navsldmixer,'UserData');
up1=get(handles.navsldbrowse,'Position');
up2=get(handles.navsldmixer,'Position');

up1(2)=up(2)-1;
up1(4)=round(ud1(4)/(ud1(4)+ud2(4))*up(4));
up2(2)=up1(2)+up1(4);
up2(4)=up(4)-up1(4)+2;
set(handles.navsldbrowse,'Position',up1);
set(handles.navsldmixer,'Position',up2);

% exit button is tied to right edge
ud=get(handles.navbtexit,'UserData');
up=ud;
up(2)=ud(2)+DY;
set(handles.navbtexit,'Position',up);


% --- Executes on button press in navrbtshadow.
function navrbtshadow_Callback(hObject, eventdata, handles)
% hObject    handle to navrbtshadow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of navrbtshadow
global data
data.shadow=get(hObject,'Value');

% only need to redraw image with new selection
utdgetimg(handles);
utdshow(handles);


% --------------------------------------------------------------------
function uimshadowmbr_Callback(hObject, eventdata, handles)
% hObject    handle to uimshadowmbr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.uimshadowmbr,'Checked'),'off'))
    set(handles.uimshadowmbr,'Checked','on');
else
    set(handles.uimshadowmbr,'Checked','off');
end

% only need to redraw image with new selection
utdgetimg(handles);
utdshow(handles);


% --------------------------------------------------------------------
function uimmult_Callback(hObject, eventdata, handles)
% hObject    handle to uimmult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.uimmult,'Checked'),'off'))
    set(handles.uimmult,'Checked','on');
else
    set(handles.uimmult,'Checked','off');
end


% only need to redraw image with new selection
utdgetimg(handles);
utdshow(handles);


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cat data proof

% free up memory
data.overhead=7e7;

% smoothing factor
smooth=10;

% downsize factor
factor=10;


imap=get(handles.navrbtmap,'Value');
id=data.ids;

set(handles.navtxtinfo,'String','generating thumbnail...'); 
set(handles.navtxtinfo,'fontweight','normal');   
pause(0.1);

thumb=imresize(cat{1},1/factor);
thumb=zeros([size(thumb),length(cat)],'uint16');

id=proof.pmap(id+1);
for k=1:length(cat)    
    set(handles.navtxtinfo,'String',sprintf('generating... reading %i...',k));
    pause(0.01);
    tmp=proof.pmap(ifget('cat',k)+1);  tmp(~ismember(tmp,id))=0;
    [x,y]=find(tmp>0);
    box={max(1,min(x)-factor):min(size(tmp,1),max(x)+factor),...
        max(1,min(y)-factor):min(size(tmp,2),max(y)+factor)};
    
    tmp(box{:})=imdilate(tmp(box{:}),ones(factor));
    thumb(:,:,k)=imresize(tmp,1/factor,'nearest');
end
thumb=thumb(end:-1:1,:,:);

if(isfield(data,'figthumb') && ~isempty(data.figthumb))
    try 
        x=get(data.figthumb,'Position')
        figure(data.figthumb); 
    catch
        data.figthumb=figure; 
        x=get(data.figthumb,'Position'); x(3:4)=[450,350];
        set(data.figthumb,'Position',x,'menubar','none',...
            'toolbar','figure','name','thumbnail3D')
    end;
else
    data.figthumb=figure;
    x=get(data.figthumb,'Position'); x(3:4)=[450,350];
    set(data.figthumb,'Position',x,'menubar','none',...
        'toolbar','figure','name','thumbnail3D')
end


set(handles.navtxtinfo,'String','drawing thumbnail...'); 
set(handles.navtxtinfo,'fontweight','normal');   
pause(0.1);

data.thumb=patch(isosurface(thumb>0,0.95));
reducepatch(data.thumb,0.1);
data.thumbcaps=patch(isocaps(thumb>0,0.95));
reducepatch(data.thumbcaps,0.1);


set(handles.navtxtinfo,'String','coloring thumbnail...'); 
set(handles.navtxtinfo,'fontweight','normal');   
pause(0.1);
data.thumblight=[camlight('right'),camlight('left')];
lighting gouraud
material dull

[X,Y,Z]=meshgrid(1:size(thumb,2),1:size(thumb,1),1:size(thumb,3));
R=im2double(reshape(data.cmap(thumb(:)+1,1),size(thumb)));
G=im2double(reshape(data.cmap(thumb(:)+1,2),size(thumb)));
B=im2double(reshape(data.cmap(thumb(:)+1,3),size(thumb)));
isocolors(X,Y,Z,R,G,B,data.thumb);
isocolors(X,Y,Z,R,G,B,data.thumbcaps);
set(data.thumb,'FaceColor','interp','EdgeColor','none');
set(data.thumbcaps,'FaceColor','interp','EdgeColor','none');

xlabel('x'); ylabel('y'); zlabel('height'); title('thumbnail3D');
set(handles.navtxtinfo,'String','done'); 
set(handles.navtxtinfo,'fontweight','normal');   
pause(0.1);
    
data.overhead=3e7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wc=ifget(var,ik)
% this function implements interface to slices in file-by-file format
wc=[];

global al cat debug proof data

 % check that we don't have it already
if(strcmp(var,'al') & ~isempty(al{ik})) wc=al{ik}; return; end
if(strcmp(var,'cat') & ~isempty(cat{ik})) wc=cat{ik}; return; end
if(strcmp(var,'wcat') & length(debug.ind1)>=ik ...
        & ~isempty(debug.ind1{ik})) wc=debug.ind1{ik}; return; end
if(strcmp(var,'alcat') & ~isempty(al{ik}) & ~isempty(cat{ik})) return; end

% identify file name
prefix=debug.prefix;
sname=sprintf('%s.%.3i.mat',prefix,ik);

if(~exist(sname))
    fprintf('DATA FILE for section %i could not be found!!!\n',ik);
    return
end


% make sure debug is OK
debug.ind1{length(cat)+1}=[];

% build list of nonempty entries in al
indk=[];
for k=1:length(cat)
    if(~isempty(al{k}) | ~isempty(cat{k}) | ~isempty(debug.ind1{k}))
        indk=[indk,k];
    end;
end


% access only images
if(strcmp(var,'al'))
    flg=1;
    while(flg)
        try                     % try to allocate memory
            overhead=zeros(1,data.overhead);
            wc_test=zeros(data.shape,'uint8');
            clear wc_test wc1_test overhead
            flg=0;
        catch                   % if error, need to clear some memory
            if(isempty(indk))
                fprintf('WARNING: unable to allocate memory buffer!\n');
                flg=0;
%                 return
            else
                % clear most remote entry
                tmp=abs(indk-ik);
                [tmp,i]=max(tmp);

                al{indk(i)}=[];
                cat{indk(i)}=[];
                debug.ind1{indk(i)}=[];
                indk=setdiff(indk,indk(i));
            end
        end
    end

    wc=load(sname,'al'); wc=wc.al;
    al{ik}=wc;
end

if(strcmp(var,'cat'))               % access segmentation slices
    flg=1;
    while(flg)
        try                     % try to allocate memory
            overhead=zeros(1,data.overhead);
            wc1_test=zeros(data.shape,'uint32');
            clear wc_test wc1_test overhead
            flg=0;
        catch                   % if error, need to clear some memory
            if(isempty(indk))
                fprintf('WARNING: unable to allocate memory buffer!\n');
                flg=0;
%                 return
            else
                % clear most remote entry
                tmp=abs(indk-ik);
                [tmp,i]=max(tmp);

                al{indk(i)}=[];
                cat{indk(i)}=[];
                debug.ind1{indk(i)}=[];
                indk=setdiff(indk,indk(i));
            end
        end
    end
        
    wc=load(sname,'cat'); wc=wc.cat;

    % apply draw-editing alterations
    if(length(data.tmp{20})>=ik)
        if(~isempty(data.tmp{20}{ik}))        
            wc(data.tmp{20}{ik}{1})=data.tmp{20}{ik}{2};
        end
    end

    % data-size reduction for cat
    if(data.mmax==1) wc=logical(wc);
    elseif(data.mmax<2^8) wc=uint8(wc);
    elseif(data.mmax<2^16) wc=uint16(wc); end
    
    cat{ik}=wc;
end

if(strcmp(var,'wcat'))               % access segmentation slices
    flg=1;
    while(flg)
        try                     % try to allocate memory
            overhead=zeros(1,data.overhead);
            wc_test=zeros(data.shape,'uint32');
            clear wc_test wc1_test overhead
            flg=0;
        catch                   % if error, need to clear some memory
            if(isempty(indk))
                fprintf('WARNING: unable to allocate memory buffer!\n');
                flg=0;
%                 return
            else
                % clear most remote entry
                tmp=abs(indk-ik);
                [tmp,i]=max(tmp);

                al{indk(i)}=[];
                cat{indk(i)}=[];
                debug.ind1{indk(i)}=[];
                indk=setdiff(indk,indk(i));
            end
        end
    end
        
    wc=load(sname,'wcat'); wc=wc.wcat;

    % data-size reduction for cat
    if(data.mmax==1) wc=logical(wc);
    elseif(data.mmax<2^8) wc=uint8(wc);
    elseif(data.mmax<2^16) wc=uint16(wc);  end
    
    debug.ind1{ik}=wc;
end
   

if(strcmp(var,'alcat'))                % access images
    flg=1;
    while(flg)
        try                     % try to allocate memory
            overhead=zeros(1,data.overhead);
            wc_test=zeros(data.shape,'uint8');
            wc1_test=zeros(data.shape,'uint32');
            clear wc_test wc1_test overhead
            flg=0;
        catch                   % if error, need to clear some memory
            if(isempty(indk))
                fprintf('WARNING: unable to allocate memory buffer!\n');
                flg=0;
%                 return
            else
                % clear most remote entry
                tmp=abs(indk-ik);
                [tmp,i]=max(tmp);

                al{indk(i)}=[];
                cat{indk(i)}=[];
                indk=setdiff(indk,indk(i));
            end
        end
    end

    wc=load(sname,'al','cat'); wc1=wc.cat; wc=wc.al;

    % apply draw-editing transformation
    if(length(data.tmp{20})>=ik)
        if(~isempty(data.tmp{20}{ik}))
            wc1(data.tmp{20}{ik}{1})=data.tmp{20}{ik}{2};
        end
    end

    % data-size reduction for cat
    if(data.mmax==1) wc1=logical(wc1);
    elseif(data.mmax<2^8) wc1=uint8(wc1);
    elseif(data.mmax<2^16) wc1=uint16(wc1); end

    al{ik}=wc;
    cat{ik}=wc1;
    
    wc=[];
end
