%% 注：无SDR 采用ARQ
clc;
clear;
close all;

%% 参数改变
MCSnum = 4;
fullnum = 500;


%% 设置MSDU和MPDU参数
msduLength = 2304; % MSDU长度
bitsPerOctet = 8;
msdubitsLen = msduLength * bitsPerOctet; % MSDU比特流的长度 
data = round(rand(msdubitsLen*fullnum,1));
numMSDUs = ceil((length(data))/msdubitsLen); % MSDU的数目
padZeros = msdubitsLen-mod(length(data)+8,msdubitsLen);
txData = [data;zeros(padZeros,1)]; % 补零

%% 生成FCS校验位
generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

%% 生成MPDU帧首部
ACK_MACHeader =   ['D4';'00'];
MACHeader =   ['08';'02';
               '00';'00';
               'FF';'FF';'FF';'FF';'FF';'FF';
               '00';'12';'34';'56';'78';'9B';
               '00';'12';'34';'56';'78';'9B';
               '00';'00'];
lengthMACheader = 24; 
lengthACKMACheader = 2;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
sequenceindex = 23;
txDataBits = zeros(0,1);

%% 设置802.11a基带信号
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.MCS = MCSnum;                   % 调制方式
nonHTcfg.ChannelBandwidth = 'CBW20';
nonHTcfg.NumTransmitAntennas = 1;   % 天线数量
nonHTcfg.PSDULength = lengthMPDU;   % PSDU长度

%% 设置接收机参数
msduLength = 2304;
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.MCS = MCSnum;                   % Modulation: 64QAM Rate: 2/3
nonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna
nonHTcfg.ChannelBandwidth = 'CBW20';
chanBW = nonHTcfg.ChannelBandwidth;

nonHTcfg.PSDULength = msduLength+28;   % Set the PSDU length
bitsPerOctet = 8;
lengthMACheader = 24; 
lengthACKMACheader = 2;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
lengthACKMPDU = lengthACKMACheader+lengthFCS;

ACKnonHTcfg = wlanNonHTConfig;         % Create packet configuration
ACKnonHTcfg.MCS = MCSnum;                   % Modulation: 64QAM Rate: 2/3
ACKnonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna
ACKnonHTcfg.ChannelBandwidth = 'CBW20';
ACKnonHTcfg.PSDULength = lengthACKMPDU;

sequenceindex = 23;
fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % Transmit sample rate in MHz
osf = 1.5;
ACKFlag = 0;
seqnum = 0;

%% 得到PSDU数据
indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

number = 0;
errrecord = [];

for NoiseSNR = 8:0.2:16
errnum = 0;
%% 发送与接收过程
for ind = 1:numMSDUs
    %% 成MPDU
    frameBody = txData((ind-1)*msdubitsLen+1:msdubitsLen*ind,:); % 生成MSDU
    frameHeader = MACHeader;
    frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
    FCS = fcsGenerator([frameHeaderBits;frameBody]);
    frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
    txDataBits = [frameHeaderBits;frameBody;frameFCS];

    %% 生成基带波形
    txWaveform = wlanWaveformGenerator(txDataBits,nonHTcfg);

    fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 基带采样率
    osf = 1.5;                     % 过采样率
    BasebandSampleRate = fs.*osf;
    txWaveform  = resample(txWaveform,fs*osf,fs); % 升采样
    fprintf('\n生成第%d个WLAN发射数据\n',ind);
    txWaveform1 = [txWaveform;zeros(10000,1)];
    %% 模拟发射端到接收端信道
    txWaveform = awgn(txWaveform1,NoiseSNR);

    %% 设置接收端接收信号长度
    captureLength_RX = length(txWaveform);

    %% 接收端开始接收数据
    burstCaptures_RX = txWaveform(1:captureLength_RX,1);
    
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

    %% 接收端处理数据
    while (searchOffset + minPktLen) <= rxWaveformLen    
        pktOffset = wlanPacketDetect(rxWaveform, chanBW, searchOffset, 0.9); %检测到包后的位置偏移

        pktOffset = searchOffset+pktOffset; % 调整偏移量
        if isempty(pktOffset) || (pktOffset+double(indLSIG(2))>rxWaveformLen)
            errnum = errnum + 1;
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
        fprintf('\n发射端第%d个包在序列%d被检测到\n',ind,pktOffset+1);

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

        refSym = wlanClosestReferenceSymbol(eqSym,rxNonHTcfg);

        % 生成MPDU的FCS
        generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
        fcsGenerator = comm.CRCGenerator(generatorPolynomial);
        fcsGenerator.InitialConditions = 1;
        fcsGenerator.DirectMethod = true;
        fcsGenerator.FinalXOR = 1;

        % 将PSDU进行分解
        rxPSDU = double(rxPSDU);
        RXMACHeader = rxPSDU(1:lengthMACheader*bitsPerOctet,1); % 得到包的头部
        rxBit = rxPSDU(lengthMACheader*bitsPerOctet+1:...
                end-lengthFCS*bitsPerOctet,1); % 得到包的数据
        RXMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1); % 得到包的FCS校验位
        sequence = rxPSDU((sequenceindex-1)*bitsPerOctet+1:sequenceindex*bitsPerOctet,1);
        packetSeq = bi2de(sequence'); % 得到包的序号;    
        rxBitMatrix(:,ind) = rxBit;

        searchOffset = pktOffset+double(indNonHTData(2)); % 完成一个包的检测，更改偏移量
        % FCS检测
        FCS = fcsGenerator([RXMACHeader;rxBit]);
        TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
        if(RXMACFCS == TrueFCS)
            break;
        else
            errnum = errnum+1;
            break;
        end
    end
end

%% 接收端解码过程
% rxData = rxBitMatrix(:);
% rxData = rxData(1:length(data));
% bitErrorRate = comm.ErrorRate;
% number = number+1;
% err = bitErrorRate(rxData,txData(1:length(rxData)));
% errrecord(number) = err(1);
% 
number = number+1;
PER(number) = errnum/fullnum;
SNRrecord(number) = NoiseSNR;
if(errnum == 0)
    break;
end

end
plot(SNRrecord,PER)