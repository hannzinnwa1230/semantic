function RTSFlag = RTSACKDetect(burstCaptures_TX)
%%
% function：receive and detect RTS ACK signal
% input：
% burstCaptures_TX：received signal
% output: 
% RTSFlag：RTS receive flag

global fcsGenerator;
%% preparation
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.ChannelBandwidth = 'CBW20';
chanBW = nonHTcfg.ChannelBandwidth;

indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

RTSFlag = 0;

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth);
osf = 1.5;

RTSref = [0;0;1;0;1;1;0;1];
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

%% process data
while (searchOffset + minPktLen) <= ACKrxWaveformLen    
    pktOffset = wlanPacketDetect(ACKrxWaveform, chanBW, searchOffset, 0.9); 

    pktOffset = searchOffset+pktOffset; 
    if (isempty(pktOffset) || pktOffset+double(indLSIG(2))>ACKrxWaveformLen)
        disp('** TX not detect RTS ACK **');
        break;
    end

    % 进行粗频率偏移矫正和符号同步
    nonHT = ACKrxWaveform(pktOffset+(indLSTF(1):indLSIG(2))); 
    coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); 
    nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); 
    fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW);
    pktOffset = pktOffset+fineTimingOffset; 

    if (pktOffset<0) || ((pktOffset+minPktLen)>ACKrxWaveformLen) 
        searchOffset = pktOffset+1.5*lstfLen; 
        continue; 
    end
    fprintf('\nTX RST ACK detected in %d\n',pktOffset+1);

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
    disp(['  ' format ' format chek success']);
    if ~strcmp(format,'Non-HT')
        fprintf('  a not Non-HT format detected\n');
        searchOffset = pktOffset+1.5*lstfLen;
        continue;
    end

    [recLSIGBits,failCheck] = wlanLSIGRecover(...
           nonHT(indLSIG(1):indLSIG(2),:),chanEstLLTF,noiseVarNonHT,chanBW);
    if failCheck
        fprintf('  L-SIG check fail \n');
        searchOffset = pktOffset+1.5*lstfLen;
        continue; 
    else
        fprintf('  L-SIG check sucess \n');
    end

    [lsigMCS,lsigLen,rxSamples] = helperInterpretLSIG(recLSIGBits,sr);

    if (rxSamples+pktOffset)>length(ACKrxWaveform)
        disp('** not enough sample to check **');
        break;
    end

    ACKrxWaveform(pktOffset+(1:rxSamples),:) = helperFrequencyOffset(...
        ACKrxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

    rxNonHTcfg = wlanNonHTConfig;
    rxNonHTcfg.MCS = lsigMCS;
    rxNonHTcfg.PSDULength = lsigLen;
    rxNonHTcfg.ChannelBandwidth = nonHTcfg.ChannelBandwidth;

    indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data'); 

    [rxPSDU,~] = wlanNonHTDataRecover(...
           ACKrxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), ...
           chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

    rxPSDU = double(rxPSDU);
    ACKMACHeader = rxPSDU(1:lengthACKMACheader*bitsPerOctet,1); 
    ACKMACcontent = rxPSDU(lengthACKMACheader*bitsPerOctet+1:end-lengthFCS*bitsPerOctet);
    ACKMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1); 

    searchOffset = pktOffset+double(indNonHTData(2)); 

    FCS = fcsGenerator([ACKMACHeader;ACKMACcontent]);
    TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
    if(ACKMACFCS == TrueFCS)
        if(isequal(ACKMACcontent(1:8),RTSref))
            disp('  RTS ACK CRC check success');
            RTSFlag = 1;
        end
        break;
    else
        disp('  RTS ACK CRC check fail');
        continue;
    end
end  


end