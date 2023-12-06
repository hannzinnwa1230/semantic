function DataReceive(burstCaptures_RX,windowlength)
%%
% function：receive and detect data
% input：
% burstCaptures_TX：received signal
% windowlength: the length of window

%% global variable
global SNRsum;
global rxwindowind;
global rxBitwindow;
global rxBitcount;
global rxBitind
global window;
global seqnum;

global constellation;
global fcsGenerator;
global evmCalculator;

%% preparation
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.ChannelBandwidth = 'CBW20';
chanBW = nonHTcfg.ChannelBandwidth;

indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth);
osf = 1.5;

lengthFCS = 4; 
lengthMACheader = 24;
bitsPerOctet = 8;
sequenceindex = 23;

rxWaveform = resample(burstCaptures_RX,fs,fs*osf); % 降采样
rxWaveformLen = size(rxWaveform,1);
searchOffset = 0; 

lstfLen = double(indLSTF(2)); 
minPktLen = lstfLen;

sr = fs; 
fineTimingOffset = []; 
packetSeq = [];
SNRsum = 0;

%% process data
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
    fprintf('\n接收端第%d个窗口在序列%d被检测到\n',rxBitind+1,pktOffset+1);

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
end