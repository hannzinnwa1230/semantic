function  MCSnum = MCSmode2num(MCSmode)
%% 
% 功能：根据MCSmode决定MCSnum
% input：
% MCSmode
% output：
% MCSnum
   
if(MCSmode == 0)
    MCSnum = [6,1];
else
    MCSnum = [5,4];
end

end