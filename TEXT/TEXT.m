%% %%%%% SDR发送端 %%%%%
clc;
clear;
close all;

%% 申明全局变量
global fcsGenerator;
global ACKFailnum;  
global numMSDUs;
global MCSnum;

% global ACKFlag;
% global ACKFailnum;    
% global ARQFlag;
% global ARQcount;

%% 参数改变
MCSnum = [6,1]; % 调制方式
alpha = 25; % 水印嵌入因子
windowlength = 5; % 滑动窗口大小
fc = 3.232e9; % 载波频率
fc_ack = 1.732e9;


%% 编码处理
img = imread("lena.tiff");
% img = imresize(img,[256,256]);
watermark = imread("uestc.bmp");
[img_watermark,realwatermark] = AddWatermark(img,watermark,alpha);
data = DataForm(img_watermark);

%% 设置MSDU和MPDU参数
msduLength = 2304; % MSDU长度
bitsPerOctet = 8;
msdubitsLen = msduLength * bitsPerOctet; % MSDU比特流的长度 
numMSDUs = ceil((length(data)+8)/msdubitsLen); % MSDU的数目
num_bi = de2bi(numMSDUs,8)';
padZeros = msdubitsLen-mod(length(data)+8,msdubitsLen);
txData = [num_bi;data;zeros(padZeros,1)]; % 补零
windownum = ceil(numMSDUs/windowlength);

%% 生成FCS校验位
generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

%% 生成MPDU帧首部
RTS_MACHeader = ['B4';'00';'00';'00']; % 00101101
lengthMACheader = 24; 
lengthRTSMACheader = 4;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
lengthRTSMPDU = lengthRTSMACheader+lengthFCS;

%% 设置802.11a基带信号
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.ChannelBandwidth = 'CBW20';
nonHTcfg.NumTransmitAntennas = 1;   % 天线数量
nonHTcfg.PSDULength = lengthMPDU;   % PSDU长度
chanBW = nonHTcfg.ChannelBandwidth;

RTSnonHTcfg = wlanNonHTConfig;         % Create packet configuration
RTSnonHTcfg.ChannelBandwidth = 'CBW20';
RTSnonHTcfg.NumTransmitAntennas = 1;   % 天线数量
RTSnonHTcfg.PSDULength = lengthRTSMPDU;   % PSDU长度

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 基带采样率
osf = 1.5;                     % 过采样率
BasebandSampleRate = fs.*osf;

%% 初始化SDR
deviceNameSDR = 'Pluto';
radio = sdrdev(deviceNameSDR); 
txGain = -10;
powerScaleFactor = 0.8; % 发射因子

%% 生成RTS信号数据
RTSframeHeader = RTS_MACHeader;
RTSframeHeaderBits = reshape((de2bi(hex2dec(RTSframeHeader)))',[],1);
FCS = fcsGenerator(RTSframeHeaderBits);
RTSframeFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
RTSBits = [RTSframeHeaderBits;RTSframeFCS];
RTSnonHTcfg.MCS = 0;
RTSWaveform = wlanWaveformGenerator(RTSBits,RTSnonHTcfg);
RTSWaveform  = resample(RTSWaveform,fs*osf,fs);

%% 循环发射RTS信号并等待回应
RTSFlag = 0;
while(RTSFlag == 0)
    %% 设置SDR发射机
    sdrTransmitter = sdrtx(deviceNameSDR); % Transmitter properties
    sdrTransmitter.RadioID = 'usb:0';

    sdrTransmitter.BasebandSampleRate = BasebandSampleRate; 
    sdrTransmitter.CenterFrequency = fc; % 2.5G频段
    sdrTransmitter.ShowAdvancedProperties = true;
    sdrTransmitter.Gain = txGain;
    
    %% 发射RTS信号
    fprintf('\n发射RTS信号\n');
    RTSWaveform = RTSWaveform.*(1/max(abs(RTSWaveform))*powerScaleFactor); % 防止发射机饱和
    sdrTransmitter.transmitRepeat(RTSWaveform);

    pause(2)
    release(sdrTransmitter);

    %% 设置发射机为接收状态
    sdrTransmitter = sdrrx('Pluto');
    sdrTransmitter.RadioID = 'usb:0';

    sdrTransmitter.BasebandSampleRate = BasebandSampleRate;
    sdrTransmitter.CenterFrequency = fc_ack;
    sdrTransmitter.GainSource = 'AGC Slow Attack';
    sdrTransmitter.OutputDataType = 'double';       
    
    %% 发射端开始接收数据
    captureLength_TX = 100000;
    burstCaptures_TX = capture(sdrTransmitter, captureLength_TX, 'Samples');
    RTSFlag = RTSACKDetect(burstCaptures_TX);
end

%% 发射数据
for ind = 1:windownum
    ACKFailnum = 0;
    ARQFlag = 1;
    ARQcount = [];
    while(ARQFlag == 1)
        if(ACKFailnum == 2)
            MCSnum = [5,4];
        end
        if(isequal(MCSnum,[6,1]))
             fprintf('发射机调制方式为64QAM\n'); 
        else
             fprintf('发射机调制方式为16QAM\n'); 
        end
        %% 数据处理
        txWaveform = WlanGenerate(txData,ind,ARQcount,windowlength);
        fprintf('\n生成第%d个WLAN发射数据\n',ind);  

        %% 设置SDR发射机
        sdrTransmitter = sdrtx(deviceNameSDR); % Transmitter properties
        sdrTransmitter.RadioID = 'usb:0';

        sdrTransmitter.BasebandSampleRate = BasebandSampleRate; 
        sdrTransmitter.CenterFrequency = fc; % 2.5G频段
        sdrTransmitter.ShowAdvancedProperties = true;
        sdrTransmitter.Gain = txGain;

        %% 发射信号
        txWaveform = txWaveform.*(1/max(abs(txWaveform))*powerScaleFactor); % 防止发射机饱和
        sdrTransmitter.transmitRepeat(txWaveform);
%         pause(5);
        pause(20);
        release(sdrTransmitter)

        %% 设置发射机为接收状态
        sdrTransmitter = sdrrx('Pluto');
        sdrTransmitter.RadioID = 'usb:0';

        sdrTransmitter.BasebandSampleRate = BasebandSampleRate;
        sdrTransmitter.CenterFrequency = fc_ack;
        sdrTransmitter.GainSource = 'AGC Slow Attack';
        sdrTransmitter.OutputDataType = 'double';       

        %% 设置发射端接收信号长度
        captureLength_TX = 1000000;
        ACKFlag = 0;
        while(ACKFlag == 0)
            %% 接收端开始接收数据
            burstCaptures_TX = capture(sdrTransmitter, captureLength_TX, 'Samples');
           [ACKFlag,ARQFlag,ARQcount] = ACKDetect(burstCaptures_TX,ind);
        end
    end
end