function varargout = pinpu_main(varargin)
% PINPU_MAIN MATLAB code for pinpu_main.fig
%      PINPU_MAIN, by itself, creates a new PINPU_MAIN or raises the existing
%      singleton*.
%
%      H = PINPU_MAIN returns the handle to a new PINPU_MAIN or the handle to
%      the existing singleton*.
%
%      PINPU_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PINPU_MAIN.M with the given input arguments.
%
%      PINPU_MAIN('Property','Value',...) creates a new PINPU_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pinpu_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pinpu_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pinpu_main

% Last Modified by GUIDE v2.5 08-Jan-2024 19:34:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pinpu_main_OpeningFcn, ...
                   'gui_OutputFcn',  @pinpu_main_OutputFcn, ...
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


% --- Executes just before pinpu_main is made visible.
function pinpu_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pinpu_main (see VARARGIN)

    
% Choose default command line output for pinpu_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pinpu_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pinpu_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in input_signal.
function input_signal_Callback(hObject, eventdata, handles)
% hObject    handle to input_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
global ecg;
global fs;
% signal=load('ecg.txt');%导入心电文件
%  t = signal(:, 1);%第一列是时间
% ecg = signal(:, 2);%第二列是幅值
fs = 360;%采样频率360Hz
axes(handles.axes1) % handles.axes1即为坐标轴axes1的句柄
plot(t,ecg);
axis([0,3.5,-0.8,1.5]);
title("心电信号时域示意图");xlabel('时间s')
xlabel('时间/s');ylabel('幅值/mV');

% --- Executes on button press in yuanpinpu_show.
function yuanpinpu_show_Callback(hObject, eventdata, handles)
% hObject    handle to yuanpinpu_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
global fs;
global ecg;
global ww;
ECG=fft(ecg);
L=length(t);
ww=(0:L-1)/L*fs;
axes(handles.axes2);%绘制坐标系
plot(ww,abs(ECG),'linewidth',1);
axis([0,180,0,100]);
title("心电信号频域示意图");xlabel('频率Hz')
xlabel('频率/Hz');ylabel('幅值/db');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in lvboxuanze.
function lvboxuanze_Callback(hObject, eventdata, handles)
% hObject    handle to lvboxuanze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lvboxuanze contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lvboxuanze

var=get(handles.lvboxuanze,'value');%可以更改value值

axes(handles.axes1);%绘制坐标系
switch var
case 1
        global ecg
        global t
        global ww;
        global f1
        Hz = str2double(get(handles.edit3, 'String')); 
        %带阻滤波器滤除工频干扰        
        fp=[Hz-1 Hz+1];fs=[Hz-0.5 Hz+0.5];
%         fp=[59 61];fs=[59.5 60.5];
        Rp = 3;        % 通带最大衰减（dB）
        As = 18;       % 阻带最大衰减（dB）
        Fs=360;
        wap=fp/(Fs/2);%通带截止频率，范围在0~1之间，1对应采样频率的一半
        was=fs/(Fs/2);%阻带截止频率
        [n,wn]=buttord(wap,was,Rp,As);
        [bz,az]=butter(n,wn,'stop');
        [h1,w]=freqz(bz,az,3600,Fs);
        f1 = filter(bz, az, ecg);%滤波!滤除工频干扰
        %显示滤波后的时域波形
        axes(handles.axes1) % handles.axes1即为坐标轴axes1的句柄
        plot(t,f1,'linewidth',1);
        axis([0,3.5,-0.8,1.5]);
        title('滤除工频干扰后时域图');
        xlabel('时间/s');ylabel('幅值/mV');
        %显示滤波后的频域波形
        y1=fft(f1);
        axes(handles.axes2) % handles.axes1即为坐标轴axes1的句柄
        plot(ww,abs(y1),'linewidth',1);
        title('滤除工频干扰后频谱图');
        xlabel('频率/Hz');ylabel('幅值/db');
        axis([0,180,0,100])
case 2
        global ecg
        global t
        global ww;
        global f1
        global f2
        %%利用高通滤除基线漂移
        Hz = str2double(get(handles.edit3, 'String')); 
        fp=Hz+1;fs=Hz;
%         fp=2;fs=1;
        Rp=3;As=18;
        Fs=360;
        wap=fp/(Fs/2);
        was=fs/(Fs/2);
        [N,Wn]=buttord(wap,was,Rp,As );
        [bz2,az2]=butter(N,Wn,'high');
        [h2,w2]=freqz(bz2,az2,3600,Fs);
        f2 = filter(bz2, az2, f1);%滤除基线漂移
        %显示滤波后的时域波形
        axes(handles.axes1) % handles.axes1即为坐标轴axes1的句柄
        plot(t,f2,'linewidth',1);
        title('滤除基线漂移后时域图');
        axis([0,3.5,-0.8,1.5]);
        xlabel('时间/s');ylabel('幅值/mV');
        %显示滤波后的频域波形
        y1=fft(f2);
        axes(handles.axes2) % handles.axes1即为坐标轴axes1的句柄
        plot(ww,abs(y1),'linewidth',1);
        title('滤除基线漂移后频谱图');
        xlabel('频率/Hz');ylabel('幅值/db');
        axis([0,180,0,100])
        ecg = f2;
end



% --- Executes during object creation, after setting all properties.
function lvboxuanze_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lvboxuanze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in choose_1.
function choose_1_Callback(hObject, eventdata, handles)
% hObject    handle to choose_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
global ecg
FileName= uigetfile({'*.txt'},'选择文件');%导入文件
str = FileName;
set(handles.edit2,'string',str);
data = load(str);
t = data(:, 1);
ecg = data(:, 2);


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in gonglvpu.
function gonglvpu_Callback(hObject, eventdata, handles)
% hObject    handle to gonglvpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gonglvpu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gonglvpu
var=get(handles.gonglvpu,'value');%可以更改value值
global t
global ecg
global fs
axes(handles.axes3);%绘制坐标系
nfft = 512;   % fft计算点数
switch var
case 1
    % 周期图法
    window=boxcar(length(ecg)); %矩形窗
    [psd1,f]=periodogram(ecg,window,nfft,fs); %直接法
    psd1 = psd1 / max(psd1);
    plot(f,10*log10(psd1+0.000001));title('周期图法');grid on
    xlabel('频率(Hz)'); % 添加x轴标签
    ylabel('功率谱密度(dB/Hz)'); % 添加y轴标签
    %归一化通常通过除以PSD的最大值来实现，
    %这样可以使得PSD的最大值等于1。
    %有时候，为了保持数值的稳定性，会在分母上加上一个很小的常数（例如0.000001），以避免除以零的情况。

case 2
    % 自相关法
    cxn=xcorr(ecg,'unbiased'); %计算序列的自相关函数
    CXk=fft(cxn,nfft);
    psd2=abs(CXk);
    index=0:round(nfft/2-1);
    k=index*fs/nfft;%计算对应的频率值k
    psd2 = psd2 / max(psd2);
    psd2=10*log10(psd2(index+1)+0.000001);
    plot(k,psd2);title('自相关法');grid on
    xlabel('频率(Hz)'); % 添加x轴标签
    ylabel('功率谱密度(dB/Hz)'); % 添加y轴标签
case 3
% 计算信号长度
sigLength = length(ecg);
% 设置FFT计算点数
nfft = 512; % FFT点数，应该是2的幂
% 设置窗函数和重叠点数
window = hamming(sigLength); % 使用汉明窗
noverlap = 0.75*sigLength; % 重叠点数，通常是窗长的50%到75%
% 使用pwelch函数进行功率谱估计
[Pxx, f] = pwelch(ecg, window, noverlap, nfft, fs);
% 绘制功率谱密度
plot(f, 10*log10(Pxx));
xlabel('频率 (Hz)');
ylabel('功率谱密度 (dB/Hz)');
title('心电信号的Welch功率谱密度估计');
grid on;
case 4
    % 自回归（AR模型）
    psd3 = pyulear(ecg, fs, nfft); %使用pyulear函数计算心电信号的AR谱估计。
    psd3=psd3/max(psd3);
    psd3 = psd3 / max(psd3);
    index=0:round(nfft/2-1);
    k=index*fs/nfft;%计算对应的频率值k
    psd3=10*log10(psd3(index+1)+0.000001);
    plot(k, psd3);title('AR谱估计');grid on;
    xlabel('频率(Hz)'); % 添加x轴标签
    ylabel('功率谱密度(dB/Hz)'); % 添加y轴标签
end

% --- Executes during object creation, after setting all properties.
function gonglvpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gonglvpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 小波变换
% 导入ECG信号数据
global t
global ecg
% 设置采样率
fs = 360; % 采样率，单位：Hz
% 连续小波变换
wavename='cmor3-3';
totalscal=50;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(ecg,scals,wavename); % 求连续小波系数
axes(handles.axes4);%绘制坐标系
imagesc(t,f,abs(coefs));
set(gca,'YDir','normal')
colorbar;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
title('小波时频图');



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
