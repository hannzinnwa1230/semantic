function databits = DataForm(image)
%%
% 功能：将输入图片压缩形成发射数据比特流
% input：
% image：输入图片
% output：
% databit：输出数据比特流
% coderate：压缩比

%% 进行图片压缩
data = JPEGEncode(image);
% data = RLCEncode(datain);

%% huffman编码
len = length(data);
uniquedata = unique(data);
for i = 1:length(uniquedata)
    num(i,1) = length(find(data == uniquedata(i,1)));
    p(1,i) = num(i,1)/len;
end
dict = huffmandict(uniquedata,p);
datacode = huffmanenco(data,dict);

%% 记录长度
len = reshape(de2bi(length(datacode),32),[],1);
databits = datacode;

end