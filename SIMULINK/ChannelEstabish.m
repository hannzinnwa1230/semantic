function signalout = ChannelEstabish(signalin,SNR)
%%
% 功能：建立多径信道模型
% input：
% signalin：输入信号
% SNR：信噪比
% output: 
% signalin：输出信号

fs = 20e6; % Hz
pathDelays = [0 1 5 50 100]*1e-9; % sec
avgPathGains = [0 -5.9 -10 -12.3 -15.6]; % dB
fD = 2; % Hz
rayleighchan = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays, ...
    'AveragePathGains',avgPathGains, ...
    'MaximumDopplerShift',fD);
% tgnChan = wlanTGnChannel('SampleRate',fs, ...
%     'LargeScaleFadingEffect','shadowing');
% signalout = tgnChan(signalin);
signalout = awgn(rayleighchan(signalin),SNR);

end