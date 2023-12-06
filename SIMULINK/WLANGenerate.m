function [dataout,BasebandSampleRate] = WLANGenerate(data,MCSnum)
%%
% 功能：将输入数据经过WLAN调制
% input：
% data：输入数据
% output：
% dataout：输出数据

%% 设置MSDU和MPDU参数
msduLength = 2304; % MSDU长度
bitsPerOctet = 8;
msdubitsLen = msduLength * bitsPerOctet; % MSDU比特流的长度 
numMSDUs = ceil((length(data)+8)/msdubitsLen); % MSDU的数目
num_bi = de2bi(numMSDUs,8)';
padZeros = msdubitsLen-mod(length(data)+8,msdubitsLen);
txData = [num_bi;data;zeros(padZeros,1)]; % 补零
save txData;

%% 生成FCS校验位
generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
fcsGenerator = comm.CRCGenerator(generatorPolynomial);
fcsGenerator.InitialConditions = 1;
fcsGenerator.DirectMethod = true;
fcsGenerator.FinalXOR = 1;

%% 生成MPDU帧首部
MACHeader =   ['08';'02';
               '00';'00';
               'FF';'FF';'FF';'FF';'FF';'FF';
               '00';'12';'34';'56';'78';'9B';
               '00';'12';'34';'56';'78';'9B';
               '00';'00'];
lengthMACheader = 24; 
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
sequenceindex = 23;
txDataBits = zeros(0,1);

%% 生成MSDU和MPDU
for ind = 0:numMSDUs-1 
    frameBody = txData(ind*msdubitsLen+1:msdubitsLen*(ind+1),:); % 生成MSDU
    Sequence = dec2hex(ind,2);
    frameHeader = MACHeader;
    frameHeader(sequenceindex,:) = num2str(Sequence);
    frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
    FCS = fcsGenerator([frameHeaderBits;frameBody]);
    frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); % 取CRC校验码后32位
    txDataBits = [txDataBits;[frameHeaderBits;frameBody;frameFCS]];
end

%% 设置802.11a基带信号
nonHTcfg = wlanNonHTConfig;         % Create packet configuration
nonHTcfg.MCS = MCSnum;                   % 调制方式
nonHTcfg.ChannelBandwidth = 'CBW20';
nonHTcfg.NumTransmitAntennas = 1;   % 天线数量
nonHTcfg.PSDULength = lengthMPDU;   % PSDU长度
scramblerInitialization = randi([1 127],numMSDUs,1); % 设置扰码

%% 生成基带波形
txWaveform = wlanWaveformGenerator(txDataBits,nonHTcfg, ...
    'NumPackets',numMSDUs,'IdleTime',20e-6, ...
    'ScramblerInitialization',scramblerInitialization);

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 基带采样率
osf = 1.5;                     % 过采样率
BasebandSampleRate = fs.*osf;
txWaveform  = resample(txWaveform,fs*osf,fs); % 升采样
dataout = txWaveform;
fprintf('\n生成WLAN发射数据\n')

end