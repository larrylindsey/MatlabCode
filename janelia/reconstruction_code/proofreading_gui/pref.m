function varargout = pref(varargin)
% PREF M-file for pref.fig
%      PREF, by itself, creates a new PREF or raises the existing
%      singleton*.
%
%      H = PREF returns the handle to a new PREF or the handle to
%      the existing singleton*.
%
%      PREF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREF.M with the given input arguments.
%
%      PREF('Property','Value',...) creates a new PREF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pref_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pref_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pref

% Last Modified by GUIDE v2.5 13-Jul-2006 17:59:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pref_OpeningFcn, ...
                   'gui_OutputFcn',  @pref_OutputFcn, ...
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


% --- Executes just before pref is made visible.
function pref_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pref (see VARARGIN)

% Choose default command line output for pref
handles.output = hObject;

% use data to exchange preferences info
global data

% assign rbt-group
set(handles.rbtcompactmode,'Value',data.pref.light);
set(handles.rbtclickselection,'Value',data.pref.selmode);
set(handles.rbtshowdebug,'Value',data.pref.dbg);
set(handles.rbtmajorselect,'Value',data.pref.selmajor);
set(handles.rbtundoclean,'Value',data.pref.undocleanup);
set(handles.rbtsinglewatershed,'Value',data.pref.watershed);
set(handles.rbtdefconn,'Value',data.pref.con==4);
set(handles.rbtcleanbuf,'Value',data.pref.delbuf);
set(handles.rbtdouble,'value',data.pref.ccdouble);

% set edt-group
str=get(handles.ppcolormap,'string');
l=false(1,length(str));
if(~isstr(data.pref.clmp)) l(end)=1; 
else
    for i=1:length(l) l(i)=strcmp(data.pref.clmp,str{i}); end
end
if(max(l)==0) l(end)=1; end
i=find(l,1);
set(handles.ppcolormap,'Value',i);

set(handles.edtpensize,'string',num2str(data.pref.edrange));
set(handles.edtshading,'string',...
    [num2str(data.pref.shade(1)),',',num2str(data.pref.shade(2))]);
set(handles.edtmovein,'string',num2str(data.pref.w));
set(handles.edttodosize,'string',num2str(data.pref.todolength));
set(handles.ppdftpensize,'value',data.pref.dftpen);
set(handles.ppdftcorrected,'value',data.pref.dftcor+1);
set(handles.ppdftquicklink,'value',data.pref.dftcklink+1);
str=sprintf(',%i',data.pref.objdskt); str=str(2:end);
set(handles.edtthrmode,'string',str);

str='';
for i=1:length(data.pref.altsel)
    if(strcmp(data.pref.altsel{i},'normal')) str=[str,',left']; end
    if(strcmp(data.pref.altsel{i},'extend')) str=[str,',shift']; end
    if(strcmp(data.pref.altsel{i},'alt')) str=[str,',right(ctrl)']; end
end
str=str(2:end);
strX=get(handles.ppselection,'string');
l=false(1,6);
for i=1:6 l(i)=strcmp(str,strX{i}); end
i=find(l,1);
set(handles.ppselection,'value',i);

% hotkeys selections
str=get(handles.pphotkeys,'string');

% list of hot-keys labels
% Xstr=cell(1,length(str));
% for i=1:length(str)
%     ind=findstr('.',str{i});
%     Xstr{i}=str{i}(1:ind(1)-1);
% end
Xstr=data.pref.keydscrp;

% reodering indexes
rids=1:length(data.pref.keys);
A=data.pref.keys;
A=A(rids);

% list of hotkeys
Ystr=cell(1,length(A));
for i=1:length(A)
    if(A(i)==30) Ystr{i}='up';
    elseif(A(i)==31) Ystr{i}='down';
    elseif(A(i)==28) Ystr{i}='left';
    elseif(A(i)==29) Ystr{i}='right';
    else Ystr{i}=char(A(i)); end
end

% form output string
for i=1:length(A)
    a=length(Xstr{i})+length(Ystr{i});
    b=max(0,14-a);
    s=repmat('.',[1,b]);
    str{i}=[Xstr{i},s,Ystr{i}];
end
set(handles.pphotkeys,'string',str);
    
    
    
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pref wait for user response (see UIRESUME)
% uiwait(handles.prefmain);


% --- Outputs from this function are returned to the command line.
function varargout = pref_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edtmovein_Callback(hObject, eventdata, handles)
% verify input for move in size box
x=str2num(get(hObject,'string'));
if(length(x)~=1)
    set(handles.txtinfo,'string','invalid move-in size numerical');
    set(handles.txtinfo,'fontweight','bold');
    set(hObject,'String','');
else
    set(handles.txtinfo,'string','GUI preferences editor');
    set(handles.txtinfo,'fontweight','normal');
end

% --- Executes during object creation, after setting all properties.
function edtmovein_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtmovein (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtshading_Callback(hObject, eventdata, handles)
% verify input at shading
x=str2num(get(hObject,'string'));
if(length(x)==0)
    set(handles.txtinfo,'string','invalid shading numericals');
    set(handles.txtinfo,'fontweight','bold');
    set(hObject,'String','');
else
    set(handles.txtinfo,'string','GUI preferences editor');
    set(handles.txtinfo,'fontweight','normal');
end


% --- Executes during object creation, after setting all properties.
function edtshading_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtshading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppcolormap.
function ppcolormap_Callback(hObject, eventdata, handles)
% hObject    handle to ppcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppcolormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppcolormap


% --- Executes during object creation, after setting all properties.
function ppcolormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppdftpensize.
function ppdftpensize_Callback(hObject, eventdata, handles)
% hObject    handle to ppdftpensize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppdftpensize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppdftpensize


% --- Executes during object creation, after setting all properties.
function ppdftpensize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppdftpensize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in rbtclickselection.
function rbtclickselection_Callback(hObject, eventdata, handles)
% hObject    handle to rbtclickselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtclickselection


% --- Executes on button press in rbtshowdebug.
function rbtshowdebug_Callback(hObject, eventdata, handles)
% hObject    handle to rbtshowdebug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtshowdebug



function edtthrmode_Callback(hObject, eventdata, handles)
% verify thr-mode setup input
x=str2num(get(hObject,'string'));
if(length(x)==0)
    set(handles.txtinfo,'string','invalid thr mode setup numerical');
    set(handles.txtinfo,'fontweight','bold');
    set(hObject,'String','');
else
    set(handles.txtinfo,'string','GUI preferences editor');
    set(handles.txtinfo,'fontweight','normal');
end



% --- Executes during object creation, after setting all properties.
function edtthrmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtthrmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pphotkeys.
function pphotkeys_Callback(hObject, eventdata, handles)
% reassings hotkey
set(handles.txtinfo,'string','press key associated with this function');
set(handles.txtinfo,'fontweight','bold');
set(handles.pphotkeys,'userdata',1);
% move focus
uicontrol(handles.btquit);

function keydown(hObject,eventdata,handles)
% catch key pressed

if(isempty(get(handles.pphotkeys,'userdata'))) return; end

% what was clicked
s=double(get(handles.prefmain,'CurrentCharacter'));
if(isempty(s)) return; end

str=get(handles.pphotkeys,'string');
ipos=get(handles.pphotkeys,'value');

Xstr=str{ipos};
k=findstr('.',Xstr);
Xstr=Xstr(1:k-1);


if(s==30) Ystr='up';
elseif(s==31) Ystr='down';
elseif(s==28) Ystr='left';
elseif(s==29) Ystr='right';
else Ystr=char(s); end

% form output string
a=length(Xstr)+length(Ystr);
b=max(0,15-a);
s=repmat('.',[1,b]);
str{ipos}=[Xstr,s,Ystr];

set(handles.pphotkeys,'string',str);
set(handles.txtinfo,'string','GUI preferences editor');
set(handles.txtinfo,'fontweight','normal');
set(handles.pphotkeys,'userdata',[]);




% --- Executes during object creation, after setting all properties.
function pphotkeys_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pphotkeys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbtsinglewatershed.
function rbtsinglewatershed_Callback(hObject, eventdata, handles)
% hObject    handle to rbtsinglewatershed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtsinglewatershed


% --- Executes on button press in rbtundoclean.
function rbtundoclean_Callback(hObject, eventdata, handles)
% hObject    handle to rbtundoclean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtundoclean


% --- Executes on button press in rbtmajorselect.
function rbtmajorselect_Callback(hObject, eventdata, handles)
% hObject    handle to rbtmajorselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtmajorselect



function edttodosize_Callback(hObject, eventdata, handles)
% verify todo input
x=str2num(get(hObject,'string'));
if(length(x)~=1)
    set(handles.txtinfo,'string','invalid todo page length numerical');
    set(handles.txtinfo,'fontweight','bold');
    set(hObject,'String','');
else
    set(handles.txtinfo,'string','GUI preferences editor');
    set(handles.txtinfo,'fontweight','normal');
end


% --- Executes during object creation, after setting all properties.
function edttodosize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edttodosize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppdftcorrected.
function ppdftcorrected_Callback(hObject, eventdata, handles)
% hObject    handle to ppdftcorrected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppdftcorrected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppdftcorrected


% --- Executes during object creation, after setting all properties.
function ppdftcorrected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppdftcorrected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% verify pen-size input
x=str2num(get(hObject,'String'));
if(length(x)~=1)
    set(handles.txtinfo,'string','invalid pen-size numerical');
    set(handles.txtinfo,'fontweight','bold');
    set(hObject,'String','');
else
    set(handles.txtinfo,'string','GUI preferences editor');
    set(handles.txtinfo,'fontweight','normal');
end


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in rbtcleanbuf.
function rbtcleanbuf_Callback(hObject, eventdata, handles)
% hObject    handle to rbtcleanbuf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtcleanbuf



function edtselection_Callback(hObject, eventdata, handles)
% verify click-selection sequence
str=get(hObject,'string');
flg=~isempty(strfind('left',str)) & ~isempty(strfind('shift',str)) ...
    & (~isempty(strfind('right',str)) | ~isempty(strfind('ctrl',str)));

if(~flg)
    set(handles.txtinfo,'string','invalid click-sequence');
    set(handles.txtinfo,'fontweight','bold');
    set(hObject,'String','');
else
    set(handles.txtinfo,'string','GUI preferences editor');
    set(handles.txtinfo,'fontweight','normal');
end



% --- Executes during object creation, after setting all properties.
function edtselection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btquit.
function btquit_Callback(hObject, eventdata, handles)
% implements export of preferences

% closereq

% use data to exchange preferences info
global data

% process rbt-group
data.pref.light=get(handles.rbtcompactmode,'Value');
data.pref.selmode=get(handles.rbtclickselection,'Value');
data.pref.dbg=get(handles.rbtshowdebug,'Value');
data.pref.selmajor=get(handles.rbtmajorselect,'Value');
data.pref.undocleanup=get(handles.rbtundoclean,'Value');
data.pref.watershed=get(handles.rbtsinglewatershed,'Value');
a=get(handles.rbtdefconn,'Value');
if(a) data.pref.con=4; else data.pref.con=8; end
data.pref.delbuf=get(handles.rbtcleanbuf,'Value');
data.pref.ccdouble=get(handles.rbtdouble,'value');

% process edt-group
str=get(handles.ppcolormap,'string'); iPos=get(handles.ppcolormap,'value');
str=str{iPos};
if(~strcmp('other',str)) 
    if(~strcmp(data.pref.clmp,str)) data.cmap=[]; end
    data.pref.clmp=str; 
end


str=get(handles.edtpensize,'string');
a=str2num(str);
if(~isempty(a)) data.pref.edrange=a(1); end

str=get(handles.edtshading,'string');
a=str2num(str);
if(length(a)>1) data.pref.shade=a(1:2); 
    elseif(length(a)==1) data.pref.shade=a*[1 0.6]; end

str=get(handles.edtmovein,'string');
a=str2num(str);
if(~isempty(a)) data.pref.w=a(1); end

str=get(handles.edttodosize,'string');
a=str2num(str);
if(~isempty(a)) data.pref.todolength=a(1); end

data.pref.dftpen=get(handles.ppdftpensize,'value');

data.pref.dftcor=get(handles.ppdftcorrected,'value')-1;
data.pref.dftcklink=get(handles.ppdftquicklink,'value')-1;

str=get(handles.edtthrmode,'string');
a=str2num(str);
if(~isempty(a) | isempty(str)) data.pref.objdskt=a; end

str=get(handles.ppselection,'string');
i=get(handles.ppselection,'value');
str=str{i};
k=1;
while(~isempty(str))
    i=findstr(',',str);
    if(~isempty(i)) 
        i=i(1); strX=str(1:i-1); 
    else
        strX=str; i=length(str); 
    end
    if(strcmp(strX,'left')) data.pref.altsel{k}='normal'; end
    if(strcmp(strX,'shift')) data.pref.altsel{k}='extend'; end
    if(~isempty(strfind(strX,'right')) | ~isempty(strfind(strX,'ctrl'))) 
        data.pref.altsel{k}='alt'; 
    end
    str=str(i+1:end);
    k=k+1;
end


% hotkeys selections
str=get(handles.pphotkeys,'string');

% list of hot-keys labels
Xstr=cell(1,length(str));
for i=1:length(str)
    ind=findstr('.',str{i});
    Xstr{i}=str{i}(1:ind(1)-1);
    Ystr{i}=str{i}(ind(end)+1:end);
end

for i=1:length(Ystr)
    if(strcmp(Ystr{i},'up')) A(i)=30;
    elseif(strcmp(Ystr{i},'down')) A(i)=31;
    elseif(strcmp(Ystr{i},'left')) A(i)=28;
    elseif(strcmp(Ystr{i},'right')) A(i)=29;
    else A(i)=uint8(Ystr{i}); end
end

% reodering indexes
rids=1:length(str);
data.pref.keys=A(rids);

prefs=data.pref;
save pref.mat prefs

closereq;

data.pref.selmode=~data.pref.selmode;
gui('navrbtselstyle_Callback',data.fmain,[],guidata(data.fmain))




% --- Executes on selection change in ppselsequence.
function ppselsequence_Callback(hObject, eventdata, handles)
% hObject    handle to ppselsequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppselsequence contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppselsequence


% --- Executes during object creation, after setting all properties.
function ppselsequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppselsequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in rbtdouble.
function rbtdouble_Callback(hObject, eventdata, handles)
% hObject    handle to rbtdouble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtdouble


