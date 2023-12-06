function RTSFlag = RTSDetect(burstCaptures_RX)
%%
% function：receive and detect RTS signal
% input：
% burstCaptures_TX：received signal
% output: 
% RTSFlag：RTS receive flag

global fcsGenerator;
%% preparation
nonHTcfg = wlanNonHTConfig;         
nonHTcfg.ChannelBandwidth = 'CBW20';
chanBW = nonHTcfg.ChannelBandwidth;

indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); 
Ns = indLSIG(2)-indLSIG(1)+1; 

RTSFlag = 0;

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth);
osf = 1.5;

RTSref = [0;0;1;0;1;1;0;1];
lengthRTSMACheader = 4;
lengthFCS = 4; 
bitsPerOctet = 8;

rxWaveform = resample(burstCaptures_RX,fs,fs*osf); 
rxWaveformLen = size(rxWaveform,1);
searchOffset = 0; 

lstfLen = double(indLSTF(2)); 
minPktLen = lstfLen;

sr = fs; 
fineTimingOffset = []; 

%% process data
while (searchOffset + minPktLen) <= rxWaveformLen    
    pktOffset = wlanPacketDetect(rxWaveform, chanBW, searchOffset, 0.9); 

    pktOffset = searchOffset+pktOffset;
    if isempty(pktOffset) || (pktOffset+double(indLSIG(2))>rxWaveformLen)
        disp('** 接收端没有检测到RTS包 **');
        break;
    end

    nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:); 
    coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW);
    nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset);
    fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); 
    pktOffset = pktOffset+fineTimingOffset; 

    if (pktOffset<0) || ((pktOffset+minPktLen)>rxWaveformLen) 
        searchOffset = pktOffset+1.5*lstfLen; 
        continue; 
    end
    fprintf('\n接收端RTS信号在序列%d被检测到\n',pktOffset+1);

    nonHT = rxWaveform(pktOffset+(1:7*Ns),:); 
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
    if (rxSamples+pktOffset)>length(rxWaveform)
        disp('** 没有足够的样点去检测包 **');
        break;
    end

    rxWaveform(pktOffset+(1:rxSamples),:) = helperFrequencyOffset(...
        rxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

    rxNonHTcfg = wlanNonHTConfig;
    rxNonHTcfg.MCS = lsigMCS;
    rxNonHTcfg.PSDULength = lsigLen;
    rxNonHTcfg.ChannelBandwidth = nonHTcfg.ChannelBandwidth;

    indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data'); 

    [rxPSDU,eqSym] = wlanNonHTDataRecover(...
           rxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), ...
           chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

    rxPSDU = double(rxPSDU);
    RTSMACHeader = rxPSDU(1:lengthRTSMACheader*bitsPerOctet,1); 
    RTSMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1); 

    searchOffset = pktOffset+double(indNonHTData(2)); 

    FCS = fcsGenerator(RTSMACHeader);
    TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
    if(RTSMACFCS == TrueFCS)
        if(isequal(RTSMACHeader(1:8),RTSref))
            disp('  RTS帧CRC检测成功');
            RTSFlag = 1;
            break;
        end
    else
        disp('  RTS帧CRC检测失败');
        continue;
    end                       
end
end