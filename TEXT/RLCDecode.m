function datacode = RLCDecode(data)
%%
% 功能：对输入数据进行行程解码
% input：
% data：输入数据(列向量)
% output：
% code：输出解码

%% 进行行程解码
[len,~] = size(data);
num = 1;
for i = 1:2:len
    L = data(i+1,1);
    for j = 1:L
        datacode(num,1) = data(i,1); 
        num = num+1;
    end
end

end