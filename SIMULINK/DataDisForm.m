function image = DataDisForm(datacode,originalimage)
%%
% 功能：将输入图片压缩形成发射数据比特流
% input：
% data：输入数据流
% output：
% image：输出图片

%% 提取图片尺寸
size_img = size(originalimage);
row = size_img(1);
col = size_img(2);

%% huffman解码
data = JPEGEncode(originalimage);
len = length(data);
uniquedata = unique(data);
for i = 1:length(uniquedata)
    num(i,1) = length(find(data == uniquedata(i,1)));
    p(1,i) = num(i,1)/len;
end
dict = huffmandict(uniquedata,p);

data = huffmandeco(datacode,dict);

%% 进行图片解压
% dataout = RLCDecode(data);
image = JPEGDecode(data,row,col);

end