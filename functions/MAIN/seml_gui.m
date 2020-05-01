function varargout = seml_gui(varargin)
% SEML_GUI M-file for seml_gui.fig
%      This is a GUI that gathers information required to run SEML and
%      places it into variable Sinfo with structure as follows:
%
%      Sinfo =
%           name: 'Lake Name'
%            lat: <value>
%           elev: <value>
%        sheight: <value>
%        mld_del: <value>

%            met: 'Path to met data'
%         tchain: 'Path to tchain data'
%           bath: 'Path to hypsographic data'
%            atn: 'Input Data' -OR- 'Average value'
%          LWout: 'Input Data' -OR- 'Calculated'
%          SWout: 'Input Data' -OR- 'Calculated'
%       comments: 'Comments'
%
%      SEML_GUI, by itself, creates a new SEML_GUI or raises the existing
%      singleton*.
%
%      H = SEML_GUI returns the handle to a new SEML_GUI or the handle to
%      the existing singleton*.
%
%      SEML_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEML_GUI.M with the given input arguments.
%
%      SEML_GUI('Property','Value',...) creates a new SEML_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seml_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seml_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help seml_gui

% Last Modified by GUIDE v2.5 26-Jan-2018 15:21:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @seml_gui_OpeningFcn, ...
    'gui_OutputFcn',  @seml_gui_OutputFcn, ...
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

% --- Executes just before seml_gui is made visible.
function seml_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seml_gui (see VARARGIN)

% Choose default command line output for seml_gui
handles = guidata(hObject);
handles.output = hObject;

handles.lakename = {'Create New Lake' 'ClearLake' 'Lawrence' 'Toolik Main' 'Mono2007','E5','Carioca'};
handles.lakeinfo =              { ''     39           42             68         38       68   19;  % Latitude
                                  ''     405         270             700        1946     700  268;  % Elevation
                                  ''     2            2             5            3        3   10};  % Sensor Height
%                    ''     180          180        1000}; % Fetch
handles.mld_name = {'d_T' 'grad_T' 'rho'};
handles.lakebath = {'';
    '\\scully\sally12\Kellogg\bathymetry\hypso\hypso_Wintergreen.mat';
    '\\scully\sally12\Kellogg\bathymetry\hypso\hypso_Lawrence.mat';
    '\\scully\sally14\Data\TOOLIK\bathymetry\hypso_toolik_main.mat';
    '\\scully\sally14\Data\Mono\Mono07\bathymetry\hypso_mono_sept_2007.mat';
    '\\scully.icess.ucsb.edu\Sally14\Data\toolik\bathymetry_other_lakes\E5_hypso.mat';''};

i = 1;  %Populate inital values with Bassett


set(handles.lakemenu,'string',handles.lakename)
set(handles.name,'string',handles.lakename{i})
set(handles.lat,'string',num2str(handles.lakeinfo{1,i}))
set(handles.elev,'string',num2str(handles.lakeinfo{2,i}))
set(handles.sheight,'string',num2str(handles.lakeinfo{3,i}))
set(handles.tchain,'string','')
%set(handles.mld_del,'string','0.02')
set(handles.met,'string','')
set(handles.bath,'string',handles.lakebath{i})
set(handles.spcond,'string','280')
set(handles.fsal,'string','0.6')
set(handles.dT,'string','0.2')
set(handles.drhodz,'string','0.003')

set(handles.loggertype,'string','RBR')


% try
%     lakstr = evalin('base','Sinfo.name;');
%     [c ia ib] = intersect(lakstr,handles.lakename);
%     set(handles.lakemenu,'value',ib)
%     lakemenu_Callback(handles.lakemenu, eventdata, handles);
%     metstr = evalin('base','Sinfo.met');
%     set(handles.met,'string',metstr)
%     tchstr = evalin('base','Sinfo.tchain');
%     set(handles.tchain,'string',tchstr)
%     handles.Sinfo.r = evalin('base','Sinfo.T');
%     handles.Sinfo.LWout = evalin('base','Sinfo.LWout');
%     handles.Sinfo.SWout = evalin('base','Sinfo.SWout');
%     tmfcn = ['warndlg([''Warning!  Paths to data have been pre-loaded.  Please check '...
%         'to ensure that they are correct.''],''Warning!'')'];
%     tm = timer('TimerFcn',tmfcn,'startdelay',.1);
%     start(tm);
% catch
% end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes seml_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = seml_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;

% --- Executes on selection change in lakemenu.
function lakemenu_Callback(hObject, eventdata, handles)
% hObject    handle to lakemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lakeind = get(hObject,'Value');
lakename = handles.lakename{lakeind};

switch lakename
    case handles.lakename(2:end)
        set(handles.name,'string',lakename)
        set(handles.lat,'string',num2str(handles.lakeinfo{1,lakeind}))
        set(handles.elev,'string',num2str(handles.lakeinfo{2,lakeind}))
        set(handles.sheight,'string',num2str(handles.lakeinfo{3,lakeind}))
        set(handles.bath,'string',handles.lakebath{lakeind})
    case handles.lakename(1)
        set(handles.name,'string','')
        set(handles.lat,'string','')
        set(handles.elev,'string','')
        set(handles.sheight,'string','')
        set(handles.bath,'string','')
end

% Hints: contents = get(hObject,'String') returns lakemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lakemenu


% --- Executes during object creation, after setting all properties.
function lakemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lakemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name as text
%        str2double(get(hObject,'String')) returns contents of name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lat_Callback(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat as text
%        str2double(get(hObject,'String')) returns contents of lat as a double

% --- Executes during object creation, after setting all properties.
function lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function spcond_Callback(hObject, eventdata, handles)
% hObject    handle to spcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spcond as text
%        str2double(get(hObject,'String')) returns contents of spcond as a double


% --- Executes during object creation, after setting all properties.
function spcond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fsal_Callback(hObject, eventdata, handles)
% hObject    handle to fsal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fsal as text
%        str2double(get(hObject,'String')) returns contents of fsal as a double


% --- Executes during object creation, after setting all properties.
function fsal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fsal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in loggertype.
function loggertype_Callback(hObject, eventdata, handles)
% hObject    handle to loggertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns loggertype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from loggertype
% f = get(handles.loggertype,'string');
% 
% switch f
%     case 'HOBO'
%         handles.Sinfo.drho = char(inputdlg('Density Gradient:',...
%                     'Input Required'));
% end
                


% --- Executes during object creation, after setting all properties.
function loggertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loggertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in tchainbutt.
function tchainbutt_Callback(hObject, eventdata, handles)
% hObject    handle to tchainbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename path] = uigetfile('*.mat','Select Thermistor Chain Data File');
f = fullfile(path,filename);
set(handles.tchain,'string',f)
tchain_Callback(hObject, eventdata, handles);


% function mld_type_Callback(hObject, eventdata, handles)
% % hObject    handle to mld_type (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of mld_type as text
% %        str2double(get(hObject,'String')) returns contents of mld_type as a double
% 
% mldind = get(hObject,'Value');
% 
% mld_type = handles.mld_name(mldind);
% 
% switch char(mld_type)
%     case handles.mld_name(1)
%         set(handles.mld_del,'string',num2str(0.02));
%         set(handles.mld_del_txt,'string','MLD Calc. Criterion [deg C]')
%     case handles.mld_name(2)
%         set(handles.mld_del,'string','');
%         set(handles.mld_del_txt,'string','MLD Calc. Criterion [Not Required]')
%     case handles.mld_name(3)
%         set(handles.mld_del,'string',num2str(0.2));
%         set(handles.mld_del_txt,'string','MLD Calc. Criterion [kg/m^3]')
% end
% 
% 
% % --- Executes during object creation, after setting all properties.
% function mld_type_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to mld_type (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

function tchain_Callback(hObject, eventdata, handles)
% hObject    handle to tchain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tchain as text
%        str2double(get(hObject,'String')) returns contents of tchain as a double

f = get(handles.tchain,'string');
tchvars = {'time' 'T' 'depth'}';
if length(f)>0 & sum(isspace(f)) ~= length(f)
    try
        varchk = who(tchvars{:},'-file',f);
        if numel(varchk) ~= numel(tchvars)
            warndlg({'Warning:  The required t-chain variables are not present';...
                '                 in the selected mat file'},'Missing Variables');
        end
        if f(1) ~= '\'
            warndlg(['Warning:  Selected path is specific to this computer''s '...
                'drive structure.  Please change the prefix from [Drive ' ...
                'Letter]:/ to \\computername\'],'Correct Path');
        end
    catch
        warndlg('Warning:  Invalid or corrupt *.mat file','Invalid File')
    end
else
end

% --- Executes during object creation, after setting all properties.
function tchain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tchain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function mld_del_Callback(hObject, eventdata, handles)
% % hObject    handle to mld_del (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of mld_del as text
% %        str2double(get(hObject,'String')) returns contents of mld_del as a double
% 
% % --- Executes during object creation, after setting all properties.
% function mld_del_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to mld_del (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

% --- Executes on button press in metbutt.
function metbutt_Callback(hObject, eventdata, handles)
% hObject    handle to metbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename path] = uigetfile('*.mat','Select Meteorological Data File');
f = fullfile(path,filename);
set(handles.met,'string',f)
met_Callback(hObject, eventdata, handles);

% --- Executes on button press in SEML_inp.
function SEML_inp_Callback(hObject, eventdata, handles)
% hObject    handle to SEML_inp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure

Sinfo.name      = get(handles.name,'string');
Sinfo.name(isspace(Sinfo.name)) = '_';
Sinfo.lat       = str2num(get(handles.lat,'string'));
Sinfo.elev      = str2num(get(handles.elev,'string'));
Sinfo.sheight   = str2num(get(handles.sheight,'string'));
%Sinfo.mld_del   = str2num(get(handles.mld_del,'string'));
Sinfo.spcond   = str2num(get(handles.spcond,'string'));
Sinfo.fsal   = str2num(get(handles.fsal,'string'));
Sinfo.dT   = str2num(get(handles.dT,'string'));
Sinfo.drhodz   = str2num(get(handles.drhodz,'string'));
Sinfo.met       = get(handles.met,'string');
Sinfo.tchain    = get(handles.tchain,'string');
Sinfo.bath      = get(handles.bath,'string');
Sinfo.loggertype = (get(handles.loggertype,'string'));
%Sinfo.atn       = handles.Sinfo.atn; % 
%Sinfo.drho       = handles.Sinfo.drho;
%Sinfo.LWout     = handles.Sinfo.LWout; % ACC commented on March 2015
%Sinfo.SWout     = handles.Sinfo.SWout; % ACC commented on March 2015
%Sinfo.comments  = get(handles.comments,'string'); % ACC commented on March 2015

handles.Sinfo = Sinfo;

if ~isdir([pwd '\tchains'])
    wa = questdlg([{'You do not appear to be sitting in the proper'}; {['directory.'...
        '  Would you still like to continue?']}],'Warning','No');
    switch wa
        case {'No' 'Cancel'}
            return
    end
end

seml(Sinfo);


function sheight_Callback(hObject, eventdata, handles)
% hObject    handle to sheight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sheight as text
%        str2double(get(hObject,'String')) returns contents of sheight as a double


% --- Executes during object creation, after setting all properties.
function sheight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sheight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elev_Callback(hObject, eventdata, handles)
% hObject    handle to elev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elev as text
%        str2double(get(hObject,'String')) returns contents of elev as a double


% --- Executes during object creation, after setting all properties.
function elev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function met_Callback(hObject, eventdata, handles)
% hObject    handle to met (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of met as text
%        str2double(get(hObject,'String')) returns contents of met as a double

f = get(handles.met,'string');
metvars = {'doy' 'AirT' 'RH' 'WS' 'SWin' 'LWin' 'SWout' 'LWout' 'atn'}';
if length(f)>0 && sum(isspace(f)) ~= length(f)
    try
        varchk = who(metvars{:},'-file',f);
        varmis = setdiff(metvars,varchk);
        if any(ismember(metvars(1:6),varmis))
            warndlg({'Warning:  One or more of the required meteorological';...
                'variables are not present in the selected mat file'},...
                'Missing Variables');
        else
            if ismember({'LWout'},varmis)
                handles.Sinfo.LWout = 'Calculated';
            else
                handles.Sinfo.LWout = 'Input Data';
            end
            if ismember({'SWout'},varmis)
                handles.Sinfo.SWout = 'Calculated';
            else
                handles.Sinfo.SWout = 'Input Data';
            end
            if ismember({'atn'},varmis)
                handles.Sinfo.atn = char(inputdlg('Average Attenuation Coefficient:',...
                    'Input Required'));
            else
                handles.Sinfo.atn = 'Input Data';
            end
            if f(1) ~= '\'
                warndlg(['Warning:  Selected path is specific to this computer''s '...
                    'drive structure.  Please change the prefix from [Drive ' ...
                    'Letter]:/ to \\computername\'],'Correct Path');
            end
        end
        handles.Sinfo;
    catch
        warndlg('Warning:  Invalid or corrupt *.mat file','Invalid File')
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function met_CreateFcn(hObject, eventdata, handles)
% hObject    handle to met (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bathbutt.
function bathbutt_Callback(hObject, eventdata, handles)
% hObject    handle to bathbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename path] = uigetfile('*.mat','Select Hypsographic Data File');
f = fullfile(path,filename);
set(handles.bath,'string',f)
bath_Callback(hObject, eventdata, handles)

function bath_Callback(hObject, eventdata, handles)
% hObject    handle to bath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bath as text
%        str2double(get(hObject,'String')) returns contents of bath as a double

f = get(handles.bath,'string');
hypvars = {'d' 'A'}';
if length(f)>0 & sum(isspace(f)) ~= length(f)
    try
        varchk = who(hypvars{:},'-file',f);
        if numel(varchk) ~= numel(hypvars)
            warndlg({'Warning:  The required hypsographic variables are not present';...
                '                 in the selected mat file'},'Missing Variables');
        end
        if f(1) ~= '\'
            warndlg(['Warning:  Selected path is specific to this computer''s '...
                'drive structure.  Please change the prefix from [Drive ' ...
                'Letter]:\ to \\computername\'],'Correct Path');
        end
    catch
        warndlg('Warning:  Invalid or corrupt *.mat file','Invalid File')
    end
else
end


% --- Executes during object creation, after setting all properties.
function bath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function comments_Callback(hObject, eventdata, handles)
% hObject    handle to comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comments as text
%        str2double(get(hObject,'String')) returns contents of comments as a double


% --- Executes during object creation, after setting all properties.
function comments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dT_Callback(hObject, eventdata, handles)
% hObject    handle to dT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dT as text
%        str2double(get(hObject,'String')) returns contents of dT as a double


% --- Executes during object creation, after setting all properties.
function dT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drhodz_Callback(hObject, eventdata, handles)
% hObject    handle to drhodz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drhodz as text
%        str2double(get(hObject,'String')) returns contents of drhodz as a double


% --- Executes during object creation, after setting all properties.
function drhodz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drhodz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
