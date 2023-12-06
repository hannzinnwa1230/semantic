function  MCSmode = SNRDecision(SNR)
%% 
% 功能：根据SNR决定调制方式
% input：
% SNR：信噪比
% output：
% MCSmode：调制方式
   
if(SNR > 4.9)
    MCSmode = 0;
else
    MCSmode = 15;
end

end