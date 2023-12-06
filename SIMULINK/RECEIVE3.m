%% %%%%% SDR接收端 %%%%%
clc;
clear;
close all;

%% 参数改变
alpha = 20; % 水印嵌入因子
windowlength = 4; % 滑动窗口大小
fc = 3.232e9; % 载波频率
fc_ack = 2.532e9;

%% 设置接收机参数
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

bitsPerOctet = 8;
fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % Transmit sample rate in MHz
osf = 1.5;
seqnum = 0;

txGain = 0;
powerScaleFactor = 0.8; % 发射因子

rxBitind = 0; % 接收窗口设置
rxwindowind = 0; 
rxBitwindow = [];
rxBitcount = [];
rxBitMatrix = [];
window = windowlength;

%% 得到PSDU数据
indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

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
    'Title', '接收基带信号频谱', ...
    'YLabel', 'Power spectral density', ...
    'Position', [69 376 800 450]);
spectrumScope.SampleRate = fs.*osf;

constellation = comm.ConstellationDiagram(...
    'Title', '接收基带信号星座图', ...
    'ShowReferenceConstellation', false, ...
    'Position', [878 376 460 460]);

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

    %% 设置接收端参数
    rxWaveform = resample(burstCaptures_RX,fs,fs*osf); % 降采样
    rxWaveformLen = size(rxWaveform,1);
    searchOffset = 0; % 初始检测偏移为0

    lstfLen = double(indLSTF(2)); % LSTF长度
    minPktLen = lstfLen*5; % 最小分组长度

    sr = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 采样率
    fineTimingOffset = []; 
    packetSeq = [];
    Index = 1;
    SNRsum = 0;

    %% 接收端处理数据
    while (searchOffset + minPktLen) <= rxWaveformLen    
        pktOffset = wlanPacketDetect(rxWaveform, chanBW, searchOffset, 0.9); %检测到包后的位置偏移

        pktOffset = searchOffset+pktOffset; % 调整偏移量
        if isempty(pktOffset) || (pktOffset+double(indLSIG(2))>rxWaveformLen)
            disp('** 接收端没有检测到RTS包 **');
            break;
        end

        % 进行粗频率偏移矫正和符号同步
        nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:); % 提取LSTF LLTF LSIG
        coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); % 载波偏移估计
        nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); % 对数据进行载波偏移矫正
        fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); % 符号偏移估计
        pktOffset = pktOffset+fineTimingOffset; 

        if (pktOffset<0) || ((pktOffset+minPktLen)>rxWaveformLen) 
            searchOffset = pktOffset+1.5*lstfLen; 
            continue; 
        end
        fprintf('\n接收端第%d个窗口在序列%d被检测到\n',rxBitind,pktOffset+1);

        % 重新提取数据
        nonHT = rxWaveform(pktOffset+(1:7*Ns),:); % 进行符号同步
        nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); % 进行载波同步

        % 提取LTF进行细频率偏移矫正和信道估计
        lltf = nonHT(indLLTF(1):indLLTF(2),:);           % 提取LTF
        fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW); 
        nonHT = helperFrequencyOffset(nonHT,fs,-fineFreqOffset); % 细频率矫正
        cfoCorrection = coarseFreqOffset+fineFreqOffset; % 总频率偏移

        lltf = nonHT(indLLTF(1):indLLTF(2),:);           % 频率矫正后再次提取LTF
        demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
        chanEstLLTF = wlanLLTFChannelEstimate(demodLLTF,chanBW); % 信道估计
        noiseVarNonHT = helperNoiseEstimate(demodLLTF); % 噪声估计

        format = wlanFormatDetect(nonHT(indLLTF(2)+(1:3*Ns),:), ...
            chanEstLLTF,noiseVarNonHT,chanBW);
        disp(['  ' format ' format检测成功']);
        if ~strcmp(format,'Non-HT')
            fprintf('  一个非Non-HT的format被检测到\n');
            searchOffset = pktOffset+1.5*lstfLen;
            continue;
        end

        %进行LSIG解码
        [recLSIGBits,failCheck] = wlanLSIGRecover(...
               nonHT(indLSIG(1):indLSIG(2),:),chanEstLLTF,noiseVarNonHT,chanBW);
        if failCheck
            fprintf('  L-SIG检测失败 \n');
            searchOffset = pktOffset+1.5*lstfLen;
            continue; 
        else
            fprintf('  L-SIG检测成功 \n');
        end

        [lsigMCS,lsigLen,rxSamples] = helperInterpretLSIG(recLSIGBits,sr);
        if (rxSamples+pktOffset)>length(rxWaveform)
            disp('** 没有足够的样点去检测包 **');
            break;
        end

        % 对整个数据包进行频率偏移矫正
        rxWaveform(pktOffset+(1:rxSamples),:) = helperFrequencyOffset(...
            rxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

        % 设置接收基带信号
        rxNonHTcfg = wlanNonHTConfig;
        rxNonHTcfg.MCS = lsigMCS;
        rxNonHTcfg.PSDULength = lsigLen;
        rxNonHTcfg.ChannelBandwidth = nonHTcfg.ChannelBandwidth;

        indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data'); % 得到PPDU的长度和位置

        % 恢复出PSDU
        [rxPSDU,eqSym] = wlanNonHTDataRecover(...
               rxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), ...
               chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

        % 绘制星座图
        constellation(reshape(eqSym,[],1)); 
        pause(0);
        release(constellation); 

        refSym = wlanClosestReferenceSymbol(eqSym,rxNonHTcfg);
        [evm.RMS,evm.Peak] = evmCalculator(refSym,eqSym);
        SNR = 10*log10(1/sqrt(evm.RMS/100));

        % 将PSDU进行分解
        rxPSDU = double(rxPSDU);
        RTSMACHeader = rxPSDU(1:lengthRTSMACheader*bitsPerOctet,1); % 得到包的头部
        RTSMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1); % 得到包的FCS校验位

        searchOffset = pktOffset+double(indNonHTData(2)); % 完成一个包的检测，更改偏移量

        % FCS检测
        FCS = fcsGenerator(RTSMACHeader);
        TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
        if(RXMACFCS == TrueFCS)
            if(isequal(RTSMACHeader(9:16),RTSref))
                disp('  RTS帧CRC检测成功');
                RTSFlag = 1;
                break;
            end
        else
            disp('  RTS帧CRC检测失败');
            continue;
        end                       
    end
    
    if(RTSFlag == 1)
        ACKframeHeader = ACK_MACHeader;
        ACKcontent = zeros(16,1);
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

        pause(5);            
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
    captureLength_RX = 1000000;

    %% 接收端开始接收数据
    burstCaptures_RX = capture(sdrReceiver, captureLength_RX, 'Samples');
    spectrumScope(burstCaptures_RX);

    %% 设置接收端参数
    rxWaveform = resample(burstCaptures_RX,fs,fs*osf); % 降采样
    rxWaveformLen = size(rxWaveform,1);
    searchOffset = 0; % 初始检测偏移为0

    lstfLen = double(indLSTF(2)); % LSTF长度
    minPktLen = lstfLen*5; % 最小分组长度

    sr = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 采样率
    fineTimingOffset = []; 
    packetSeq = [];
    Index = 1;
    SNRsum = 0;

    %% 接收端处理数据
    while (searchOffset + minPktLen) <= rxWaveformLen    
        pktOffset = wlanPacketDetect(rxWaveform, chanBW, searchOffset, 0.9); %检测到包后的位置偏移

        pktOffset = searchOffset+pktOffset; % 调整偏移量
        if isempty(pktOffset) || (pktOffset+double(indLSIG(2))>rxWaveformLen)
            disp('** 接收端没有检测到数据包 **');
            break;
        end

        % 进行粗频率偏移矫正和符号同步
        nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:); % 提取LSTF LLTF LSIG
        coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); % 载波偏移估计
        nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); % 对数据进行载波偏移矫正
        fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); % 符号偏移估计
        pktOffset = pktOffset+fineTimingOffset; 

        if (pktOffset<0) || ((pktOffset+minPktLen)>rxWaveformLen) 
            searchOffset = pktOffset+1.5*lstfLen; 
            continue; 
        end
        fprintf('\n接收端第%d个窗口在序列%d被检测到\n',rxBitind,pktOffset+1);

        % 重新提取数据
        nonHT = rxWaveform(pktOffset+(1:7*Ns),:); % 进行符号同步
        nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); % 进行载波同步

        % 提取LTF进行细频率偏移矫正和信道估计
        lltf = nonHT(indLLTF(1):indLLTF(2),:);           % 提取LTF
        fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW); 
        nonHT = helperFrequencyOffset(nonHT,fs,-fineFreqOffset); % 细频率矫正
        cfoCorrection = coarseFreqOffset+fineFreqOffset; % 总频率偏移

        lltf = nonHT(indLLTF(1):indLLTF(2),:);           % 频率矫正后再次提取LTF
        demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
        chanEstLLTF = wlanLLTFChannelEstimate(demodLLTF,chanBW); % 信道估计
        noiseVarNonHT = helperNoiseEstimate(demodLLTF); % 噪声估计

        format = wlanFormatDetect(nonHT(indLLTF(2)+(1:3*Ns),:), ...
            chanEstLLTF,noiseVarNonHT,chanBW);
        disp(['  ' format ' format检测成功']);
        if ~strcmp(format,'Non-HT')
            fprintf('  一个非Non-HT的format被检测到\n');
            searchOffset = pktOffset+1.5*lstfLen;
            continue;
        end

        %进行LSIG解码
        [recLSIGBits,failCheck] = wlanLSIGRecover(...
               nonHT(indLSIG(1):indLSIG(2),:),chanEstLLTF,noiseVarNonHT,chanBW);
        if failCheck
            fprintf('  L-SIG检测失败 \n');
            searchOffset = pktOffset+1.5*lstfLen;
            continue; 
        else
            fprintf('  L-SIG检测成功 \n');
        end

        [lsigMCS,lsigLen,rxSamples] = helperInterpretLSIG(recLSIGBits,sr);
        if (rxSamples+pktOffset)>length(rxWaveform)
            disp('** 没有足够的样点去检测包 **');
            break;
        end

        % 对整个数据包进行频率偏移矫正
        rxWaveform(pktOffset+(1:rxSamples),:) = helperFrequencyOffset(...
            rxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

        % 设置接收基带信号
        rxNonHTcfg = wlanNonHTConfig;
        rxNonHTcfg.MCS = lsigMCS;
        rxNonHTcfg.PSDULength = lsigLen;
        rxNonHTcfg.ChannelBandwidth = nonHTcfg.ChannelBandwidth;

        indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data'); % 得到PPDU的长度和位置

        % 恢复出PSDU
        [rxPSDU,eqSym] = wlanNonHTDataRecover(...
               rxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), ...
               chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

        % 绘制星座图
        constellation(reshape(eqSym,[],1)); 
        pause(0);
        release(constellation); 

        refSym = wlanClosestReferenceSymbol(eqSym,rxNonHTcfg);
        [evm.RMS,evm.Peak] = evmCalculator(refSym,eqSym);
        SNR = 10*log10(1/sqrt(evm.RMS/100));

        % 将PSDU进行分解
        rxPSDU = double(rxPSDU);
        RXMACHeader = rxPSDU(1:lengthMACheader*bitsPerOctet,1); % 得到包的头部
        rxBit = rxPSDU(lengthMACheader*bitsPerOctet+1:...
                end-lengthFCS*bitsPerOctet,1); % 得到包的数据
        RXMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1); % 得到包的FCS校验位
        sequence = rxPSDU((sequenceindex-1)*bitsPerOctet+1:sequenceindex*bitsPerOctet,1);
        packetSeq = bi2de(sequence'); % 得到包的序号   
        
        signalpow = sum((abs(rxWaveform(pktOffset+(1:indNonHTData(2))))).^2)/double((indNonHTData(2)));
        SNRshow = 10*log10(signalpow/noiseVarNonHT-1);

        % 显示EVM
        fprintf('  EVM peak: %0.3f%%  EVM RMS: %0.3f%%  SNR: %0.1fdB% ',evm.Peak,evm.RMS,SNRshow);
        fprintf('\n');
        searchOffset = pktOffset+double(indNonHTData(2)); % 完成一个包的检测，更改偏移量

        % 丢掉重复的包
        if(~isempty(rxBitwindow))
            if(ismember(packetSeq,rxBitwindow(1,:)))
                    continue;
            end
        end

        % FCS检测
        FCS = fcsGenerator([RXMACHeader;rxBit]);
        TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
        if(RXMACFCS == TrueFCS)
            disp('  MAC帧CRC检测成功');
            SNRsum = SNRsum+SNR; 
            rxwindowind = rxwindowind+1;
            rxBitwindow(:,rxwindowind) = [packetSeq;rxBit];
            rxBitcount(1,rxwindowind) = mod(packetSeq,windowlength);

            if (packetSeq == 0 && seqnum == 0)
                seqnum = bi2de(rxBit(1:8)');
            end

            if(rxBitind == floor(seqnum/windowlength))
                window = mod(seqnum,windowlength);
            else
                window = windowlength;
            end

            if(rxwindowind == window)
                break;
            end

            continue;
        else
            disp('  MAC帧CRC检测失败');
            continue;
        end                       
    end
    
    for i = 1:window
        windowref(1,i) = i-1;   %标准窗口
    end

    if(rxwindowind == window)
        SNR = SNRsum/window;
        MCSmode = SNRDecision(SNR);
        rxwindowind = 0; % 接收窗口设置
        rxBitMatrix(:,rxBitind*windowlength+1:rxBitind*windowlength+window) = rxBitwindow;
        rxBitwindow = [];
        rxBitcount = [];
        rxBitind = rxBitind+1;            
        window = windowlength;
        ACKcontent = [reshape(de2bi(MCSmode,8),[],1);1;1;1;1;1;1;1;1];
        lengthACKMPDU = lengthACKMACheader+lengthFCS+2;
    else
        ARQcount = reshape(de2bi(setdiff(windowref,rxBitcount)',8)',[],1);
        ACKcontent = [zeros(16,1);ARQcount];
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

    pause(15);            
    
    if(rxBitind == ceil(seqnum/windowlength) && seqnum ~= 0)
        break;
    end
end

%% 接收端解码过程
img = imread("lena.bmp");
img = imresize(img,[256,256]);
watermark = imread("uestc.bmp");
[img_watermark,realwatermark] = AddWatermark(img,watermark,alpha);
data = DataForm(img_watermark);

rxData = SequenceDepend(rxBitMatrix);
rxData = rxData(1:length(data));
tximage = img_watermark;
img_rec = DataDisForm(rxData,tximage);
watermark_pick = PickWatermark(img_rec,alpha);

%% 峰均信噪比和相关性
PSNR = PSNRCalc(img,img_watermark)
NC = NCCalc(realwatermark,watermark_pick)

%% 显示
figure(1)
imshow(img_watermark)
title("发送图像")
figure(2)
imshow(watermark_pick)
title("提取水印")