function varargout = startgui(varargin)
% STARTGUI is the loading script for Proofreading Tool GUI.
% Set global fLoadInd2=0 to disable loading of debug.ind2 variable to
% reduce loaded dataset size, or if "ind2 variable not found" error message
% is returned when loading.
%
% GUI v.7.2. by Yuriy Mishchenko Chklovskii Lab JFRC DEC 2006


% Last Modified by GUIDE v2.5 19-Jul-2006 12:55:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @startgui_OpeningFcn, ...
                   'gui_OutputFcn',  @startgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

% % THIS CHECKS IF FILE HAD BEEN PASSED AS PARAMETER
% global guiflg1978x
% if(~isempty(varargin) && isempty(guiflg1978x))
%     if(exist(varargin{1})==2)
%         fprintf('Starting with %s ...\n\n',varargin{1});
%         
%         clear global al cat debug proof data buffer
%         global al cat debug proof
%         pause(0.1);
%         warning off MATLAB:load:variableNotFound
%         load(varargin{1},'cat');
%         mmax=-Inf;
%         if(iscell(cat))
%             for k=1:length(cat) mmax=max(mmax,max(cat{k}(:))); end
%             if(mmax<2^8) for k=1:length(cat) cat{k}=uint8(cat{k});  end
%             elseif(mmax<2^16) for k=1:length(cat) cat{k}=uint16(cat{k}); end;
%             end
%         else
%             mmax=max(mmax,max(cat(:)));
%             if(mmax<2^8) cat=uint8(cat);
%             elseif(mmax<2^16) cat=uint16(cat);
%             end
%         end
%         load(varargin{1},'al');
%         warning on MATLAB:load:variableNotFound
%         load(varargin{1},'debug','proof');
%         gui;
%         clear global guiflg1978x
%         return
%     end
% else
%     guiflg1978x=1;
% end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% clear global guiflg1978x



% --- Executes just before startgui is made visible.
function startgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to startgui (see VARARGIN)

% THIS VAR DETERMINES WHETHER debug.ind2 WILL BE LOADED
global fLoadInd2

if(isempty(fLoadInd2)) fLoadInd2=1; end
    
% THIS DETERMINES IF THERE HAD BEEN INPUT FILE PASSED IN
global flg

% populate current directory list
dirs=dir('*.mat');
str={};
for i=1:length(dirs) 
    if(isempty(findstr('gui',dirs(i).name)) & ...
            isempty(findstr('pref.mat',dirs(i).name))) 
        str{end+1}=dirs(i).name; 
    end
end
set(handles.lbxdir,'string',str);

if(length(dirs)==0) set(handles.btstart,'enable','off'); end


% THIS CHECKS IF FILE HAD BEEN PASSED AS PARAMETER
handles.flg=0;
if(~isempty(varargin))
    for i=1:length(dirs)
        if(strcmp(dirs(i).name,varargin{1}))
            set(handles.lbxdir,'value',i);
            handles.flg=1;
        end
    end
end
            


% Choose default command line output for startgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes startgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = startgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if file had bee passed
if(handles.flg)
    btstart_Callback(handles.btstart,[],handles)
end

% Get default command line output from handles structure
varargout{1} = handles.output;

if(handles.flg) delete(handles.figure1); end


% --- Executes during object creation, after setting all properties.
function lbxdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbxdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
                    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lbxdir_Callback(hObject, eventdata, handles)
set(handles.btstart,'enable','on');

% --- Executes on button press in btstart.
function btstart_Callback(hObject, eventdata, handles)
% hObject    handle to btstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% loads data
global fLoadInd2

s=get(handles.lbxdir,'string');
i=get(handles.lbxdir,'value');
set(handles.btstart,'enable','off');

clear global al cat debug proof data buffer
global al cat debug proof
set(handles.txtmain,'string','Loading GUI data...');
pause(0.1);
%load(s{i},'al','cat','debug','proof');
warning off MATLAB:load:variableNotFound
load(s{i},'cat');
mmax=-Inf; 
if(iscell(cat))
    for k=1:length(cat) mmax=max(mmax,max(cat{k}(:))); end
    if(mmax<2^8) for k=1:length(cat) cat{k}=uint8(cat{k});  end
    elseif(mmax<2^16) for k=1:length(cat) cat{k}=uint16(cat{k}); end;
    end
else
    mmax=max(mmax,max(cat(:)));
    if(mmax<2^8) cat=uint8(cat);
    elseif(mmax<2^16) cat=uint16(cat);
    end
end
if(fLoadInd2) 
    ind2=[];
    load(s{i},'ind2');
    if(iscell(ind2))
        if(mmax<2^8) for k=1:length(ind2) ind2{k}=uint8(ind2{k});  end
        elseif(mmax<2^16) for k=1:length(ind2) ind2{k}=uint16(ind2{k}); end;
        end
    else
        if(mmax<2^8) ind2=uint8(ind2);
        elseif(mmax<2^16) ind2=uint16(ind2); end
    end
end
load(s{i},'al');
warning on MATLAB:load:variableNotFound    
load(s{i},'debug','proof');
if(fLoadInd2) debug.ind2=ind2; clear ind2; end

set(handles.txtmain,'string','Starting GUI now...');

gui;
[h,f]=gcbo;

% clear global guiflg1978x

delete(f);

function myclosereq

% clear global guiflg1978x
closereq

