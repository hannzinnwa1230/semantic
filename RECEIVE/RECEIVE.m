%% %%%%% SDR接收端 %%%%%
clc;
clear;
close all;

%% 申明全局变量
global SNRsum;
global rxwindowind; 
global rxBitwindow;
global rxBitcount;
global rxBitind
global window;
global seqnum;

global eqsymnum1;
global eqsym1;
global eqsymnum2;
global eqsym2;
global eqsymnum3;
global eqsym3;

global constellation;
global fcsGenerator;
global evmCalculator;

%% 参数改变
alpha = 25; % 水印嵌入因子
windowlength = 5; % 滑动窗口大小
fc = 3.232e9; % 载波频率
fc_ack = 1.732e9;

%% 设置接收机参数
RTSref = [0;0;1;0;1;1;0;1];
ACK_MACHeader = ['D4';'00';'00';'00']; % 00101011

msduLength = 2304;
lengthMACheader = 24; 
lengthACKMACheader = 4;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
lengthACKMPDU = lengthACKMACheader+lengthFCS+2;

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

txGain = -10;
powerScaleFactor = 0.8; % 发射因子

rxBitind = 0; % 接收窗口设置
rxwindowind = 0; 
rxBitwindow = [];
rxBitcount = [];
rxBitMatrix = [];
window = windowlength;

%% 生成FCS校验位
generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

%% 设置频谱仪和星座图
spectrumScope = dsp.SpectrumAnalyzer( ...
    'SpectrumType', 'Power density', ...
    'SpectralAverages', 10, ...
    'YLimits', [-130 0], ...
    'YLabel', 'Power spectral density', ...
    'Position', [69 176 800 450]);
spectrumScope.SampleRate = fs.*osf;

constellation = comm.ConstellationDiagram(...
    'ShowReferenceConstellation', false, ...
    'Position', [578 176 460 460]);

%% 设置EVM计算
evmCalculator = comm.EVM('AveragingDimensions',[1 2 3]);
evmCalculator.MaximumEVMOutputPort = true;

%% 等待RTS信号
RTSFlag = 0;
while(RTSFlag == 0)
    %% 设置SDR接收机
    sdrReceiver = sdrrx('Pluto');
    sdrReceiver.RadioID = 'usb:0';

    sdrReceiver.BasebandSampleRate = fs*osf;
    sdrReceiver.CenterFrequency = fc;
    sdrReceiver.GainSource = 'AGC Slow Attack';
    sdrReceiver.OutputDataType = 'double';  
    
    %% 设置接收端接收信号长度
    captureLength_RX = 5000;
    %% 接收端开始接收数据
    burstCaptures_RX = capture(sdrReceiver, captureLength_RX, 'Samples');
    RTSFlag = RTSDetect(burstCaptures_RX);
    if(RTSFlag == 1)
        ACKframeHeader = ACK_MACHeader;
        ACKcontent = [RTSref;zeros(8,1)];
        ACKframeHeaderBits = reshape((de2bi(hex2dec(ACKframeHeader)))',[],1);
        FCS = fcsGenerator([ACKframeHeaderBits;ACKcontent]);
        ACKframeFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
        ACKBits = [ACKframeHeaderBits;ACKcontent;ACKframeFCS];
        ACKnonHTcfg.MCS = 0;
        ACKWaveform = wlanWaveformGenerator(ACKBits,ACKnonHTcfg);
        ACKWaveform  = resample(ACKWaveform,fs*osf,fs);    

        %% 设置SDR接收机为发射状态
        sdrReceiver = sdrtx('Pluto');
        sdrReceiver.RadioID = 'usb:0';

        sdrReceiver.BasebandSampleRate = fs*osf; 
        sdrReceiver.CenterFrequency = fc_ack; % 2.5G频段
        sdrReceiver.ShowAdvancedProperties = true;
        sdrReceiver.Gain = txGain;

        %% 发射ACK信号
        powerScaleFactor = 0.8;
        ACKWaveform = ACKWaveform.*(1/max(abs(ACKWaveform))*powerScaleFactor); % 防止发射机饱和
        sdrReceiver.transmitRepeat(ACKWaveform);

        pause(3);      
        release(sdrReceiver);
    end
end

disp('接收到RTS信号，开始接收数据');

while(1)
    %% 设置SDR接收机
    sdrReceiver = sdrrx('Pluto');
    sdrReceiver.RadioID = 'usb:0';

    sdrReceiver.BasebandSampleRate = fs*osf;
    sdrReceiver.CenterFrequency = fc;
    sdrReceiver.GainSource = 'AGC Slow Attack';
    sdrReceiver.OutputDataType = 'double';  
    
    %% 设置接收端接收信号长度
    captureLength_RX = 1500000;
    
    eqsymnum1 = 1;
    eqsym1 = [];
    eqsymnum2 = 1;
    eqsym2 = [];  
    eqsymnum3 = 1;
    eqsym3 = [];  
    
    %% 接收端开始接收数据
    burstCaptures_RX = capture(sdrReceiver, captureLength_RX, 'Samples');
    spectrumScope(burstCaptures_RX);
    DataReceive(burstCaptures_RX,windowlength);
    
    windowref = zeros(1,length(window));
    for i = 1:window
        windowref(1,i) = i-1;   %标准窗口
    end

    if(rxwindowind == window)
        fprintf('\n共接受到%d个包，进行ACK确认 \n', rxwindowind);
        SNR = SNRsum/window;
        MCSmode = SNRDecision(SNR);
        rxwindowind = 0; % 接收窗口设置
        if(isempty(rxBitMatrix) == 1 || isequal(ismember(rxBitwindow(1,:),rxBitMatrix(1,:)),...
        zeros(1,window)))
            rxBitMatrix(:,rxBitind*windowlength+1:rxBitind*windowlength+window) = rxBitwindow;
            rxBitind = rxBitind+1;            
        end
        rxBitwindow = [];
        rxBitcount = [];
        window = windowlength;
        ACKcontent = [reshape(de2bi(MCSmode,8),[],1);1;1;1;1;1;1;1;1;reshape(de2bi(rxBitind,8),[],1)];
        lengthACKMPDU = lengthACKMACheader+lengthFCS+3;
    else
        fprintf('\n共接受到%d个包，请求重传 \n', rxwindowind);
        ARQcount = reshape(de2bi(setdiff(windowref,rxBitcount)',8)',[],1);
        ACKcontent = [zeros(8,1);ones(8,1);ARQcount];
        lengthACKMPDU = lengthACKMACheader+lengthFCS+2+length(ARQcount)/8;
    end
    
    windowref = [];
    ACKnonHTcfg.PSDULength = lengthACKMPDU;
    ACKframeHeader = ACK_MACHeader;
    ACKframeHeaderBits = reshape((de2bi(hex2dec(ACKframeHeader)))',[],1);
    FCS = fcsGenerator([ACKframeHeaderBits;ACKcontent]);
    ACKframeFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
    ACKBits = [ACKframeHeaderBits;ACKcontent;ACKframeFCS];
    ACKnonHTcfg.MCS = 0;
    ACKWaveform = wlanWaveformGenerator(ACKBits,ACKnonHTcfg);
    ACKWaveform  = resample(ACKWaveform,fs*osf,fs);    

    %% 设置SDR接收机为发射状态
    sdrReceiver = sdrtx('Pluto');
    sdrReceiver.RadioID = 'usb:0';

    sdrReceiver.BasebandSampleRate = fs*osf; 
    sdrReceiver.CenterFrequency = fc_ack; % 2.5G频段
    sdrReceiver.ShowAdvancedProperties = true;
    sdrReceiver.Gain = txGain;

    %% 发射ACK信号
    powerScaleFactor = 0.8;
    ACKWaveform = ACKWaveform.*(1/max(abs(ACKWaveform))*powerScaleFactor); % 防止发射机饱和
    sdrReceiver.transmitRepeat(ACKWaveform);

%     pause(5);       
    pause(30);
    release(sdrReceiver);
    
    if(rxBitind == ceil(seqnum/windowlength) && seqnum ~= 0)
        break;
    end
end

%% 接收端解码过程
img = imread("lena.tiff");
% img = imresize(img,[256,256]);
watermark = imread("uestc.bmp");
[img_watermark,realwatermark] = AddWatermark(img,watermark,alpha);
data = DataForm(img_watermark);

rxData = SequenceDepend(rxBitMatrix);
rxData = rxData(1:length(data));
tximage = img_watermark;
img_rec = DataDisForm(rxData,tximage);
watermark_pick = PickWatermark(img_rec,alpha);

%% 峰均信噪比和相关性
PSNR = PSNRCalc(img,img_watermark);
NC = NCCalc(realwatermark,watermark_pick);

%% 显示频谱和星座图
% spectrumScope(burstCaptures_RX);
% constellation(reshape(eqSym,[],1)); 
% release(constellation); 
% constellation(reshape(eqsym(:,:,3),[],1))
release(constellation); 
constellation(reshape(eqsym1(:,:,3),[],1))
%% 显示
figure(1)
subplot(221)
imshow(img_rec)
title("接收的图像")
subplot(222)
imshow(watermark_pick)
title("提取的水印")
subplot(223)
imshow(img_watermark)
title("发送的图像")
subplot(224)
imshow(realwatermark)
title("嵌入的水印")