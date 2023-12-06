function [ACKFlag,ARQFlag,ARQcount] = ACKDetect(burstCaptures_TX,ind)
%%
% function：receive and detect RTS signal
% input：
% burstCaptures_TX: received signal
% output: 
% ACKFlag：ACK receive flag
% ARQFlag: ARQ receive flag
% ARQcount:  ARQ count number
% MCSnum: MCS number

global ACKFailnum;
global MCSnum;

%% set parameter
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.ChannelBandwidth = 'CBW20';
chanBW = nonHTcfg.ChannelBandwidth;

indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth);
osf = 1.5;

lengthACKMACheader = 4;
lengthFCS = 4; 
bitsPerOctet = 8;

ACKrxWaveform = resample(burstCaptures_TX,fs,fs*osf); 
ACKrxWaveformLen = size(ACKrxWaveform,1);
searchOffset = 0; 

lstfLen = double(indLSTF(2)); 
minPktLen = lstfLen;

sr = fs; 
fineTimingOffset = []; 

ACKFlag = 0;
ARQFlag = 1;
ARQcount = [];

%% process data
while (searchOffset + minPktLen) <= ACKrxWaveformLen    
    pktOffset = wlanPacketDetect(ACKrxWaveform, chanBW, searchOffset, 0.9);

    pktOffset = searchOffset+pktOffset; 
    if (isempty(pktOffset) || pktOffset+double(indLSIG(2))>ACKrxWaveformLen)
        disp('** 发射端没有检测到ACK包 **');
        break;
    end

    nonHT = ACKrxWaveform(pktOffset+(indLSTF(1):indLSIG(2))); 
    coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); 
    nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); 
    fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); 
    pktOffset = pktOffset+fineTimingOffset; 

    if (pktOffset<0) || ((pktOffset+minPktLen)>ACKrxWaveformLen) 
        searchOffset = pktOffset+1.5*lstfLen; 
        continue; 
    end
    fprintf('\n发射端第%d个ACK包在序列%d被检测到\n',ind,pktOffset+1);

    nonHT = ACKrxWaveform(pktOffset+(1:7*Ns),:);
    nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); 

    lltf = nonHT(indLLTF(1):indLLTF(2),:);           
    fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW); 
    nonHT = helperFrequencyOffset(nonHT,fs,-fineFreqOffset); 
    cfoCorrection = coarseFreqOffset+fineFreqOffset; 

    lltf = nonHT(indLLTF(1):indLLTF(2),:);           
    demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
    chanEstLLTF = wlanLLTFChannelEstimate(demodLLTF,chanBW); 
    noiseVarNonHT = helperNoiseEstimate(demodLLTF);

    format = wlanFormatDetect(nonHT(indLLTF(2)+(1:3*Ns),:), ...
        chanEstLLTF,noiseVarNonHT,chanBW);
    disp(['  ' format ' format检测成功']);
    if ~strcmp(format,'Non-HT')
        fprintf('  一个非Non-HT的format被检测到\n');
        searchOffset = pktOffset+1.5*lstfLen;
        continue;
    end

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

    if (rxSamples+pktOffset)>length(ACKrxWaveform)
        disp('** 没有足够的样点去检测ACK包 **');
        break;
    end

    % 对整个数据包进行频率偏移矫正
    ACKrxWaveform(pktOffset+(1:rxSamples),:) = helperFrequencyOffset(...
        ACKrxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

    % 设置接收基带信号
    rxNonHTcfg = wlanNonHTConfig;
    rxNonHTcfg.MCS = lsigMCS;
    rxNonHTcfg.PSDULength = lsigLen;
    rxNonHTcfg.ChannelBandwidth = nonHTcfg.ChannelBandwidth;

    indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data'); % 得到PPDU的长度和位置

    % 恢复出PSDU
    [rxPSDU,~] = wlanNonHTDataRecover(...
           ACKrxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), ...
           chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

    % 生成MPDU的FCS
    generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
    fcsGenerator = comm.CRCGenerator(generatorPolynomial);
    fcsGenerator.InitialConditions = 1;
    fcsGenerator.DirectMethod = true;
    fcsGenerator.FinalXOR = 1;

    rxPSDU = double(rxPSDU);
    ACKMACHeader = rxPSDU(1:lengthACKMACheader*bitsPerOctet,1); 
    ACKMACcontent = rxPSDU(lengthACKMACheader*bitsPerOctet+1:end-lengthFCS*bitsPerOctet);
    ACKMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1);

    searchOffset = pktOffset+double(indNonHTData(2)); 

    FCS = fcsGenerator([ACKMACHeader;ACKMACcontent]);
    TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
    if(ACKMACFCS == TrueFCS)
        disp('  ACK帧CRC检测成功');
        ACKFlag = 1;
        if(isequal(ACKMACcontent(9:16),[1;1;1;1;1;1;1;1]) && bi2de(ACKMACcontent(17:24)') == ind)
            fprintf('\n准备传输下一个窗口 \n');
            ARQFlag = 0;
            ARQcount = [];
            MCSmode = bi2de(reshape(ACKMACcontent(1:8),1,8));
            MCSnum = MCSmode2num(MCSmode);
        elseif(isequal(ACKMACcontent(1:16),[zeros(8,1);ones(8,1)]))
            fprintf('\n准备重传 \n');
            ARQcount = bi2de(reshape(ACKMACcontent(17:end),8,[])')';
            ARQFlag = 1;
            MCSnum = [5,4];
            ACKFailnum = ACKFailnum+1;
        end
        break;
    else
        disp('  ACK帧CRC检测失败');
        ACKFlag = 0;
        ARQFlag = 1;
        continue;
    end
end
end