clc
clear
close all

%% 设置接收机参数
global numMSDUs;
global SNRsum;
global rxwindowind;
global rxBitwindow;
global rxBitcount;
global rxBitind
global window;
global seqnum;
global MCSnum;
global constellation;
global fcsGenerator;
global evmCalculator;

rxBitind = 0;
rxwindowind = 0; 
rxBitwindow = [];
rxBitcount = [];
seqnum = 0;

RTSref = [0;0;1;0;1;1;0;1];
RTS_MACHeader = ['B4';'00';'00';'00']; % 00101101
ACK_MACHeader = ['D4';'00';'00';'00']; % 00101011
MACHeader =   ['08';'02';
               '00';'00';
               'FF';'FF';'FF';'FF';'FF';'FF';
               '00';'12';'34';'56';'78';'9B';
               '00';'12';'34';'56';'78';'9B';
               '00';'00'];
msduLength = 2304;
lengthMACheader = 24; 
lengthRTSMACheader = 4;
lengthACKMACheader = 4;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
lengthACKMPDU = lengthRTSMACheader+lengthFCS+2;
lengthRTSMPDU = lengthRTSMACheader+lengthFCS;
sequenceindex = 23;

nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.ChannelBandwidth = 'CBW20';
nonHTcfg.NumTransmitAntennas = 1;   % 天线数量
nonHTcfg.PSDULength = lengthMPDU;   % PSDU长度
chanBW = nonHTcfg.ChannelBandwidth;

ACKnonHTcfg = wlanNonHTConfig;         % Create packet configuration
ACKnonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna
ACKnonHTcfg.ChannelBandwidth = 'CBW20';
ACKnonHTcfg.PSDULength = lengthACKMPDU;

bitsPerOctet = 8;
fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % Transmit sample rate in MHz
osf = 1.5;
seqnum = 0;

%% 生成FCS校验位

% %% RTS test
% ACKnonHTcfg = wlanNonHTConfig;         % Create packet configuration
% ACKnonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna
% ACKnonHTcfg.ChannelBandwidth = 'CBW20';
% ACKnonHTcfg.PSDULength = lengthACKMPDU;
% 
% ACKframeHeader = ACK_MACHeader;
% ACKcontent = [RTSref;zeros(8,1)];
% ACKframeHeaderBits = reshape((de2bi(hex2dec(ACKframeHeader)))',[],1);
% FCS = fcsGenerator([ACKframeHeaderBits;ACKcontent]);
% ACKframeFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
% ACKBits = [ACKframeHeaderBits;ACKcontent;ACKframeFCS];
% ACKnonHTcfg.MCS = 0;
% ACKWaveform = wlanWaveformGenerator(ACKBits,ACKnonHTcfg);
% ACKWaveform  = resample(ACKWaveform,fs*osf,fs);
% 
% burstCaptures_TX = [zeros(10000,1);ACKWaveform];
% RTSFlag = RTSACKDetect(burstCaptures_TX)

%% text test
MCSnum = [6,1]; % 调制方式
alpha = 45; % 水印嵌入因子
windowlength = 4; % 滑动窗口大小
window = windowlength;

%% 生成FCS校验位
generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

constellation = comm.ConstellationDiagram(...
    'Title', '接收基带信号星座图', ...
    'ShowReferenceConstellation', false, ...
    'Position', [878 376 460 460]);

%% 设置EVM计算
evmCalculator = comm.EVM('AveragingDimensions',[1 2 3]);
evmCalculator.MaximumEVMOutputPort = true;

%% 编码处理
img = imread("lena.bmp");
img = imresize(img,[256,256]);
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

ind = 1;
    ARQcount = [];
    txWaveform = WlanGenerate(txData,ind,ARQcount,windowlength);
    
rxWaveform = resample1(txWaveform,fs,fs*osf); % 降采样
% tic    
% burstCaptures_RX = [zeros(10000,1);txWaveform];
% DataReceive(burstCaptures_RX,windowlength);
% toc

% ACKcontent = [zeros(16,1);de2bi(1,8)';de2bi(3,8)';de2bi(4,8)'];
% ACKnonHTcfg.PSDULength = lengthACKMPDU+3;
% ACKframeHeader = ACK_MACHeader;
% ACKframeHeaderBits = reshape((de2bi(hex2dec(ACKframeHeader)))',[],1);
% FCS = fcsGenerator([ACKframeHeaderBits;ACKcontent]);
% ACKframeFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
% ACKBits = [ACKframeHeaderBits;ACKcontent;ACKframeFCS];
% ACKnonHTcfg.MCS = 0;
% ACKWaveform = wlanWaveformGenerator(ACKBits,ACKnonHTcfg);
% ACKWaveform  = resample(ACKWaveform,fs*osf,fs);    
% 
% ind = 1;
% burstCaptures_TX = [zeros(10000,1);ACKWaveform];
% [ACKFlag,ARQFlag,ARQcount] = ACKDetect(burstCaptures_TX,ind)

% %% 生成RTS信号数据
% RTSnonHTcfg = wlanNonHTConfig;         % Create packet configuration
% RTSnonHTcfg.ChannelBandwidth = 'CBW20';
% RTSnonHTcfg.NumTransmitAntennas = 1;   % 天线数量
% RTSnonHTcfg.PSDULength = lengthRTSMPDU;   % PSDU长度
% 
% RTSframeHeader = RTS_MACHeader;
% RTSframeHeaderBits = reshape((de2bi(hex2dec(RTSframeHeader)))',[],1);
% FCS = fcsGenerator(RTSframeHeaderBits);
% RTSframeFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
% RTSBits = [RTSframeHeaderBits;RTSframeFCS];
% RTSnonHTcfg.MCS = 0;
% RTSWaveform = wlanWaveformGenerator(RTSBits,RTSnonHTcfg);
% RTSWaveform  = resample(RTSWaveform,fs*osf,fs);
% burstCaptures_RX = [zeros(10000,1);RTSWaveform];
% 
% RTSFlag = RTSDetect(burstCaptures_RX)
