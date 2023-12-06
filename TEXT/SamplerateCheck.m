function fs = SamplerateCheck(ChannelBandwidth)
%%
% 功能：根据输入符号串判断采样率
% input：
% ChannelBandwidth：输入采样率字符串
% output: 
% fs：输出采样率

%%
switch ChannelBandwidth
    case 'CBW5'
        fs = 5e6;
    case 'CBW10'
        fs = 10e6;
    case 'CBW20'
        fs = 20e6;
    case 'CBW40'
        fs = 40e6;
    case 'CBW80'
        fs = 80e6;
    case 'CBW160'
        fs = 160e6;
end

end