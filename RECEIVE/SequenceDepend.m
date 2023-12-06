function dataout = SequenceDepend(datain)
%%
% 功能：按照包序列提取数据
% input：
% datain：输入数据
% output：
% dataout：输出数据

%% 提取顺序
Sequence = datain(1,:);

%% 按顺序提取数据
dataout = [];
maxsequence = max(Sequence);
for i = 0:maxsequence
    for j = 1:maxsequence+1
        if(datain(1,j) == i)
            if (datain(1,j) == 0)
                dataout = [dataout;datain(10:end,j)];
            else
                dataout = [dataout;datain(2:end,j)];
            end
        end
    end
end

end