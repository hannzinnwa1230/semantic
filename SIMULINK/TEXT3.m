%% %%%%% SDR发送端 %%%%%
clc;
clear;
close all;

%% 参数改变
MCSnum = [6,1]; % 调制方式
RTSMCSnum = [3,2];
alpha = 45; % 水印嵌入因子
windowlength = 4; % 滑动窗口大小
fc = 3.232e9; % 载波频率
fc_ack = 2.532e9;

%% 编码处理
img = imread("./TrainPic/lena.bmp");
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

%% 生成FCS校验位
generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

%% 生成MPDU帧首部
RTS_MACHeader = ['B4';'00';'00';'00']; % 00101101
ACK_MACHeader = ['D4';'00';'00';'00']; % 00101011
MACHeader =   ['08';'02';
               '00';'00';
               'FF';'FF';'FF';'FF';'FF';'FF';
               '00';'12';'34';'56';'78';'9B';
               '00';'12';'34';'56';'78';'9B';
               '00';'00'];
lengthMACheader = 24; 
lengthRTSMACheader = 4;
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
lengthRTSMPDU = lengthRTSMACheader+lengthFCS;
sequenceindex = 23;

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

txWaveform1 = [];
txWaveform2 = [];

%% 得到PSDU数据
indLSTF = wlanFieldIndices(nonHTcfg,'L-STF');  
indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG'); % 得到LSTF、LLTF、LSIG的长度
Ns = indLSIG(2)-indLSIG(1)+1; % OFDM符号长度

%% 初始化SDR
deviceNameSDR = 'Pluto';
radio = sdrdev(deviceNameSDR); 
txGain = 0;
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

    %% 设置发射机为接收状态
    sdrTransmitter = sdrrx('Pluto');
    sdrTransmitter.RadioID = 'usb:0';

    sdrTransmitter.BasebandSampleRate = BasebandSampleRate;
    sdrTransmitter.CenterFrequency = fc_ack;
    sdrTransmitter.GainSource = 'AGC Slow Attack';
    sdrTransmitter.OutputDataType = 'double';       

    %% 设置发射端接收长度
    captureLength_TX = 100000;
    
    %% 发射端开始接收数据
    burstCaptures_TX = capture(sdrTransmitter, captureLength_TX, 'Samples');
    
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
            disp('** 发射端没有检测到RTS的ACK包 **');
            break;
        end

        % 进行粗频率偏移矫正和符号同步
        nonHT = ACKrxWaveform(pktOffset+(indLSTF(1):indLSIG(2))); % 提取LSTF LLTF LSIG
        coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); % 载波偏移估计
        nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); % 对数据进行载波偏移矫正
        fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); % 符号偏移估计
        pktOffset = pktOffset+fineTimingOffset; 

        if (pktOffset<0) || ((pktOffset+minPktLen)>ACKrxWaveformLen) 
            searchOffset = pktOffset+1.5*lstfLen; 
            continue; 
        end
        fprintf('\n发射端第%d个ACK包在序列%d被检测到\n',ind,pktOffset+1);

        % 重新提取数据
        nonHT = ACKrxWaveform(pktOffset+(1:7*Ns),:); % 进行符号同步
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
            disp('  RTS的ACK帧CRC检测成功');
            RTSFlag = 1;
            break;
        else
            disp('  RTS的ACK帧CRC检测失败');
            continue;
        end
    end  
end
    
%% 发射数据
for ind = 1:windownum
    ACKFlag = 0;
    ACKFailnum = 0;
    ARQFlag = 1;
    ARQcount = [];
    while(ARQFlag == 1)
        if(ACKFailnum == 2)
            MCSnum = [5,4];
        end
        %% 数据处理
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
        pause(10);

        %% 设置发射机为接收状态
        sdrTransmitter = sdrrx('Pluto');
        sdrTransmitter.RadioID = 'usb:0';

        sdrTransmitter.BasebandSampleRate = BasebandSampleRate;
        sdrTransmitter.CenterFrequency = fc_ack;
        sdrTransmitter.GainSource = 'AGC Slow Attack';
        sdrTransmitter.OutputDataType = 'double';       

        %% 设置发射端接收信号长度
        captureLength_TX = 1000000;

        while(ACKFlag == 0)
            %% 接收端开始接收数据
            burstCaptures_TX = capture(sdrTransmitter, captureLength_TX, 'Samples');
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
                    disp('** 发射端没有检测到ACK包 **');
                    break;
                end

                % 进行粗频率偏移矫正和符号同步
                nonHT = ACKrxWaveform(pktOffset+(indLSTF(1):indLSIG(2))); % 提取LSTF LLTF LSIG
                coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); % 载波偏移估计
                nonHT = helperFrequencyOffset(nonHT,fs,-coarseFreqOffset); % 对数据进行载波偏移矫正
                fineTimingOffset = wlanSymbolTimingEstimate(nonHT,chanBW); % 符号偏移估计
                pktOffset = pktOffset+fineTimingOffset; 

                if (pktOffset<0) || ((pktOffset+minPktLen)>ACKrxWaveformLen) 
                    searchOffset = pktOffset+1.5*lstfLen; 
                    continue; 
                end
                fprintf('\n发射端第%d个ACK包在序列%d被检测到\n',ind,pktOffset+1);

                % 重新提取数据
                nonHT = ACKrxWaveform(pktOffset+(1:7*Ns),:); % 进行符号同步
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
                        ARQFlag = 0;
                        MCSmode = bi2de(reshape(ACKMACcontent(1:8),1,8));
                        MCSnum = MCSmode2num(MCSmode);
                    else
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
        end
    end
end