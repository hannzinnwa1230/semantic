function txWaveform = WlanGenerate(txData,ind,ARQcount,windowlength)
%%
% function：receive and detect RTS signal
% input：
% txData: text data
% ind: packet sequence
% MCSnum: MCS number 
% ARQcount：ARQ sequence
% windowlength: the length of window
% output: 
% txWaveform：text waveform

global fcsGenerator;
global numMSDUs;
global MCSnum;

%% set MSDU and MPDU parameter 
MACHeader =   ['08';'02';
               '00';'00';
               'FF';'FF';'FF';'FF';'FF';'FF';
               '00';'12';'34';'56';'78';'9B';
               '00';'12';'34';'56';'78';'9B';
               '00';'00'];
msduLength = 2304; 
lengthMACheader = 24; 
lengthFCS = 4;      
lengthMPDU = lengthMACheader+msduLength+lengthFCS;
sequenceindex = 23;
bitsPerOctet = 8;
msdubitsLen = msduLength * bitsPerOctet;

nonHTcfg = wlanNonHTConfig;        
nonHTcfg.ChannelBandwidth = 'CBW20';
nonHTcfg.NumTransmitAntennas = 1;   
nonHTcfg.PSDULength = lengthMPDU;

fs = SamplerateCheck(nonHTcfg.ChannelBandwidth); % 基带采样率
osf = 1.5;    

txWaveform1 = [];
txWaveform2 = [];

if(isempty(ARQcount) == 1)
    for windowind = 1:windowlength
        tempind = (ind-1)*windowlength+windowind;
        if(tempind <= numMSDUs)
            frameBody = txData((tempind-1)*msdubitsLen+1:msdubitsLen*tempind,:); 
            Sequence = dec2hex(tempind-1,2);
            frameHeader = MACHeader;
            frameHeader(sequenceindex,:) = num2str(Sequence);
            frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
            FCS = fcsGenerator([frameHeaderBits;frameBody]);
            frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); 
            txDataBits = [frameHeaderBits;frameBody;frameFCS];
            
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
else
    for i = 1:length(ARQcount)
        tempind = (ind-1)*windowlength+ARQcount(i)+1;
        frameBody = txData((tempind-1)*msdubitsLen+1:msdubitsLen*tempind,:);
        Sequence = dec2hex(tempind-1,2);
        frameHeader = MACHeader;
        frameHeader(sequenceindex,:) = num2str(Sequence);
        frameHeaderBits = reshape((de2bi(hex2dec(frameHeader)))',[],1);
        FCS = fcsGenerator([frameHeaderBits;frameBody]);
        frameFCS = FCS(end-lengthFCS*bitsPerOctet+1:end); 
        txDataBits = [frameHeaderBits;frameBody;frameFCS];

        nonHTcfg.MCS = MCSnum(1);
        txWaveform1 = [txWaveform1;wlanWaveformGenerator(txDataBits,nonHTcfg)];
        nonHTcfg.MCS = MCSnum(2);
        txWaveform2 = [txWaveform2;wlanWaveformGenerator(txDataBits,nonHTcfg)];
    end
end

txWaveform = [txWaveform2;txWaveform1];
txWaveform  = resample(txWaveform,fs*osf,fs); 

end