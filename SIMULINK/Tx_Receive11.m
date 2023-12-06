%% 注：无SDR 采用多编码自适应滑动窗口ARQ
clc;
clear;
close all;

%% 参数改变
MCSnum = [6,1]; % 调制方式
NoiseSNR = 25; % 信道信噪比
alpha = 25; % 水印嵌入因子
windowlength = 5; % 滑动窗口大小

%% 编码处理
img = imread("lena.tiff");
% img = imresize(img,[300,300]);
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
ACK_MACHeader = ['D4';'00';'00';'00';];
MACHeader =   ['08';'02';
               '00';'00';
               'FF';'FF';'FF';'FF';'FF';'FF';
               '00';'12';'34';'56';'78';'9B';
               '00';'12';'34';'56';'78';'9B';
               '00';'00'];
lengthMACheader = 24; 
lengthACKMACheader = 4;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
% lengthACKMPDU = lengthACKMACheader+lengthFCS;
% sequenceindex = 23;

%% 设置802.11a基带信号
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.ChannelBandwidth = 'CBW20';
nonHTcfg.NumTransmitAntennas = 1;   % 天线数量
nonHTcfg.PSDULength = lengthMPDU;   % PSDU长度

%% 设置接收机参数
msduLength = 2304;
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna
nonHTcfg.ChannelBandwidth = 'CBW20';
chanBW = nonHTcfg.ChannelBandwidth;

nonHTcfg.PSDULength = msduLength+28;   % Set the PSDU length
bitsPerOctet = 8;
lengthMACheader = 24; 
lengthACKMACheader = 4;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
lengthACKMPDU = lengthACKMACheader+1+lengthFCS;

ACKnonHTcfg = wlanNonHTConfig;         % Create packet configuration
ACKnonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna
ACKnonHTcfg.ChannelBandwidth = 'CBW20';
ACKnonHTcfg.PSDULength = lengthACKMPDU;

sequenceindex = 23;
fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % Transmit sample rate in MHz
osf = 1.5;
seqnum = 0;
txWaveform1 = [];
txWaveform2 = [];

rxBitind = 0; % 接收窗口设置
rxwindowind = 0; 
rxBitwindow = []; 
rxBitcount = [];
rxBitMatrix = [];
window = windowlength;

QAMdownFlag = 0; % QAM降低开启标志

processnum = 0;

%% 得到PSDU数据
indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

%% 发送与接收过程
for ind = 1:windownum
    ACKFlag = 0;
    ACKFailnum = 0;
    ARQFlag = 1;
    ARQcount = [];
    for windowind = 1:windowlength
        tempind = (ind-1)*windowlength+windowind;
        if(tempind <= numMSDUs)
            frameBody = txData((tempind-1)*msdubitsLen+1:msdubitsLen*tempind,:); % 生成MSDU
            Sequence = dec2hex(tempind-1,2);
            frameHeader = MACHeader;
            frameHeader(sequenceindex,:) = num2str(Sequence);
            frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
            FCS = fcsGenerator([frameHeaderBits;frameBody]);
            frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
            txDataBits = [frameHeaderBits;frameBody;frameFCS];
            %% 生成基带波形
            if(tempind == 1)
                nonHTcfg.MCS = 0;
            else
                nonHTcfg.MCS = MCSnum(1);
            end
            txWaveform1 = [txWaveform1;wlanWaveformGenerator(txDataBits,nonHTcfg)];
            nonHTcfg.MCS = MCSnum(2);
            txWaveform2 = [txWaveform2;wlanWaveformGenerator(txDataBits,nonHTcfg)];            
        else
            break;
        end
    end
    
    txWaveform = [txWaveform2;txWaveform1];
    txWaveform1 = [];
    txWaveform2 = [];
    fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 基带采样率
    osf = 1.5;                     % 过采样率
    BasebandSampleRate = fs.*osf;
    txWaveform  = resample(txWaveform,fs*osf,fs); % 升采样
    fprintf('\n生成第%d个WLAN发射数据\n',ind);
    
    while (ARQFlag == 1)
        if(isempty(ARQcount) == 1)
            for windowind = 1:windowlength
                tempind = (ind-1)*windowlength+windowind;
                if(tempind <= numMSDUs)
                    frameBody = txData((tempind-1)*msdubitsLen+1:msdubitsLen*tempind,:); % 生成MSDU
                    Sequence = dec2hex(tempind-1,2);
                    frameHeader = MACHeader;
                    frameHeader(sequenceindex,:) = num2str(Sequence);
                    frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
                    FCS = fcsGenerator([frameHeaderBits;frameBody]);
                    frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
                    txDataBits = [frameHeaderBits;frameBody;frameFCS];
                    %% 生成基带波形
                    nonHTcfg.MCS = MCSnum(1);
                    txWaveform1 = [txWaveform1;wlanWaveformGenerator(txDataBits,nonHTcfg)];
                    nonHTcfg.MCS = MCSnum(2);
                    txWaveform2 = [txWaveform2;wlanWaveformGenerator(txDataBits,nonHTcfg)];            
                else
                    break;
                end
            end
        else
            for i = 1:length(ARQcount)
                tempind = (ind-1)*windowlength+ARQcount(i)+1;
                frameBody = txData((tempind-1)*msdubitsLen+1:msdubitsLen*tempind,:); % 生成MSDU
                Sequence = dec2hex(tempind-1,2);
                frameHeader = MACHeader;
                frameHeader(sequenceindex,:) = num2str(Sequence);
                frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
                FCS = fcsGenerator([frameHeaderBits;frameBody]);
                frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
                txDataBits = [frameHeaderBits;frameBody;frameFCS];
                %% 生成基带波形
                nonHTcfg.MCS = MCSnum(1);
                txWaveform1 = [txWaveform1;wlanWaveformGenerator(txDataBits,nonHTcfg)];
                nonHTcfg.MCS = MCSnum(2);
                txWaveform2 = [txWaveform2;wlanWaveformGenerator(txDataBits,nonHTcfg)];
            end
        end
        
        txWaveform = [txWaveform2;txWaveform1];
        txWaveform1 = [];
        txWaveform2 = [];
        txWaveform  = resample(txWaveform,fs*osf,fs); % 升采样
              
        processnum = processnum + 1;
        if(processnum == 2)
            NoiseSNR = 15;
        elseif(processnum == 5)
            NoiseSNR = 25;
        end
        
        %% 模拟发射端到接收端信道
        txWaveform_s = [txWaveform;zeros(10000,1)];
        txWaveform = ChannelEstabish(txWaveform_s,NoiseSNR);

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

        %% 设置接收端接收信号长度
        captureLength_RX = length(txWaveform);

        %% 接收端开始接收数据
        burstCaptures_RX = txWaveform(1:captureLength_RX,1);
%         burstCaptures_RX = zeros(2000,1);
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
                disp('** 接收端没有检测到包 **');
                break;
            end

            % 进行粗频率偏移矫正和符号同步
            nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:); % 提取LSTF LLTF LSIG
            coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); % 载波偏移估计
            nonHT = frequencyOffset(nonHT,fs,-coarseFreqOffset); % 对数据进行载波偏移矫正
            fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); % 符号偏移估计
            pktOffset = pktOffset+fineTimingOffset; 

            if (pktOffset<0) || ((pktOffset+minPktLen)>rxWaveformLen) 
                searchOffset = pktOffset+1.5*lstfLen; 
                continue; 
            end
            fprintf('\n接收端第%d个窗口在序列%d被检测到\n',ind,pktOffset+1);

            % 重新提取数据
            nonHT = rxWaveform(pktOffset+(1:7*Ns),:); % 进行符号同步
            nonHT = frequencyOffset(nonHT,fs,-coarseFreqOffset); % 进行载波偏移矫正

            % 提取LTF进行细频率偏移矫正和信道估计
            lltf = nonHT(indLLTF(1):indLLTF(2),:);           % 提取LTF
            fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW); 
            nonHT = frequencyOffset(nonHT,fs,-fineFreqOffset); % 细频率矫正
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
            rxWaveform(pktOffset+(1:rxSamples),:) = frequencyOffset(...
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
            packetSeq = bi2de(sequence'); % 得到包的序号   
            
            signal = rxWaveform(pktOffset+(1:indNonHTData(2)));
            signalpow = sum((abs(rxWaveform(pktOffset+(1:indNonHTData(2))))).^2)/double((indNonHTData(2)));
            SNRshow = 10*log10(signalpow/noiseVarNonHT-1);
            
%             figure(2)
%             M = length(signal);
%             window = hamming(M);
%             noverlap = M/2;
%             Nfft = length(signal);
%             [Pxx2 w]=pwelch(signal,window,noverlap,Nfft);
%             f = (-length(signal)/2:length(signal)/2-1)/length(signal)*fs;
%             figure;
%             subplot(211)
%             plot(f,10*log10(Pxx2)),grid

%             signalfft = fftshift(fft(signal))/length(signal);
%             signalpsd = 10*log10(abs(signalfft).^2);
%             f = (-length(signalpsd)/2:length(signalpsd)/2-1)/length(signalpsd)*fs/1e6;
%             plot(f,signalpsd),grid
%             xlabel('频率/MHz')
%             ylabel('功率谱密度/dB')
%             title('接收基带信号功率谱')
           
            % 显示EVM
            fprintf('  EVM峰值: %0.3f%%  EVM均方根值: %0.3f%%\n',evm.Peak,evm.RMS);
            fprintf('  SNR:%0.3f%',SNRshow);
            fprintf('  \n');
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
        
               
        %% 模拟接收端到发射端信道
        ACKWaveform1 = [ACKWaveform;zeros(100000,1)];
        ACKWaveform = awgn(ACKWaveform1,NoiseSNR);
        
        %% 设置发射端接收信号长度
        captureLength_TX = length(ACKWaveform);
%         captureLength_TX = 2000;

        %% 接收端开始接收数据
        burstCaptures_TX = ACKWaveform(1:captureLength_TX,1);

       %% 设置发射端参数
        ACKrxWaveform = resample(burstCaptures_TX,fs,fs*osf); % 降采样
        ACKrxWaveformLen = size(ACKrxWaveform,1);
        searchOffset = 0; % 初始检测偏移为0

        lstfLen = double(indLSTF(2)); % LSTF长度
        minPktLen = lstfLen; % 最小分组长度

        sr = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 采样率
        fineTimingOffset = []; 
        packetSeq = [];

        %% 发射端处理数据
        while (searchOffset + minPktLen) <= ACKrxWaveformLen    
            pktOffset = wlanPacketDetect(ACKrxWaveform, chanBW, searchOffset, 0.9); %检测到包后的位置偏移

            pktOffset = searchOffset+pktOffset; % 调整偏移量
            if (isempty(pktOffset) || pktOffset+double(indLSIG(2))>ACKrxWaveformLen)
                disp('** 发射端没有检测到包 **');
                ACKFailnum = ACKFailnum+1;
                break;
            end

            % 进行粗频率偏移矫正和符号同步
            nonHT = ACKrxWaveform(pktOffset+(indLSTF(1):indLSIG(2))); % 提取LSTF LLTF LSIG
            coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); % 载波偏移估计
            nonHT = frequencyOffset(nonHT,fs,-coarseFreqOffset); % 对数据进行载波偏移矫正
            fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); % 符号偏移估计
            pktOffset = pktOffset+fineTimingOffset; 

            if (pktOffset<0) || ((pktOffset+minPktLen)>ACKrxWaveformLen) 
                searchOffset = pktOffset+1.5*lstfLen; 
                continue; 
            end
            fprintf('\n发射端第%d个ACK包在序列%d被检测到\n',ind,pktOffset+1);

            % 重新提取数据
            nonHT = ACKrxWaveform(pktOffset+(1:7*Ns),:); % 进行符号同步
            nonHT = frequencyOffset(nonHT,fs,-coarseFreqOffset); % 进行载波同步

            % 提取LTF进行细频率偏移矫正和信道估计
            lltf = nonHT(indLLTF(1):indLLTF(2),:);           % 提取LTF
            fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW); 
            nonHT = frequencyOffset(nonHT,fs,-fineFreqOffset); % 细频率矫正
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

            if (rxSamples+pktOffset)>length(ACKrxWaveform)
                disp('** 没有足够的样点去检测ACK包 **');
                break;
            end

            % 对整个数据包进行频率偏移矫正
            ACKrxWaveform(pktOffset+(1:rxSamples),:) = frequencyOffset(...
                ACKrxWaveform(pktOffset+(1:rxSamples),:),fs,-cfoCorrection);

            % 设置接收基带信号
            rxNonHTcfg = wlanNonHTConfig;
            rxNonHTcfg.MCS = lsigMCS;
            rxNonHTcfg.PSDULength = lsigLen;
            rxNonHTcfg.ChannelBandwidth = ACKnonHTcfg.ChannelBandwidth;

            indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data'); % 得到PPDU的长度和位置

            % 恢复出PSDU
            [rxPSDU,eqSym] = wlanNonHTDataRecover(...
                   ACKrxWaveform(pktOffset+(indNonHTData(1):indNonHTData(2)),:), ...
                   chanEstLLTF,noiseVarNonHT,rxNonHTcfg);

            % 生成MPDU的FCS
            generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
            fcsGenerator = comm.CRCGenerator(generatorPolynomial);
            fcsGenerator.InitialConditions = 1;
            fcsGenerator.DirectMethod = true;
            fcsGenerator.FinalXOR = 1;

            % 将PSDU进行分解
            rxPSDU = double(rxPSDU);
            ACKMACHeader = rxPSDU(1:lengthACKMACheader*bitsPerOctet,1); % 得到包的头部
            ACKMACcontent = rxPSDU(lengthACKMACheader*bitsPerOctet+1:end-lengthFCS*bitsPerOctet);
            ACKMACFCS = rxPSDU(end-lengthFCS*bitsPerOctet+1:end,1); % 得到包的FCS校验位

            searchOffset = pktOffset+double(indNonHTData(2)); % 完成一个包的检测，更改偏移量

            % FCS检测
            FCS = fcsGenerator([ACKMACHeader;ACKMACcontent]);
            TrueFCS = FCS(end-lengthFCS*bitsPerOctet+1:end);
            if(ACKMACFCS == TrueFCS)
                disp('  ACK帧CRC检测成功');
                ACKFlag = 1;
                if(isequal(ACKMACcontent(9:16),[1;1;1;1;1;1;1;1]))
                    fprintf('\n准备传输下一个窗口\n');
                    ARQFlag = 0;
                    ARQcount = [];
                    MCSmode = bi2de(reshape(ACKMACcontent(1:8),1,8));
                    MCSnum = MCSmode2num(MCSmode);
                else
                    fprintf('\n准备重传\n');
                    ARQFlag = 1;
                    ARQcount = bi2de(reshape(ACKMACcontent(17:end),8,[])')';
                    ACKFailnum = ACKFailnum+1;
                end
                break;
            else
                disp('  ACK帧CRC检测失败');
                ACKFlag = 0;
                continue;
            end
        end
        if(ACKFailnum == 2)
            MCSnum = [5,4];
%         elseif(ACKFailnum == 4)
%             nonHTcfg.MCS = 4;
%             ACKnonHTcfg.MCS = 4;
        end
    end
end

%% 接收端解码过程
rxData = SequenceDepend(rxBitMatrix);
rxData = rxData(1:length(data));

%计算误码率
% bitErrorRate = comm.ErrorRate;
% err = bitErrorRate(rxData(1:length(txData)),txData);
% fprintf('\n\n');
% fprintf('  误码率为: %0.8f%%\n',err(1));

tximage = img_watermark;
img_rec = DataDisForm(rxData,tximage);
watermark_pick = PickWatermark(img_rec,alpha);

%% 峰均信噪比和相关性
PSNR = PSNRCalc(img,img_watermark)
NC = NCCalc(realwatermark,watermark_pick)

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