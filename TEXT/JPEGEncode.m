function datacode = JPEGEncode(img)
%%
% 功能：对图片进行JPEG压缩
% input：
% img：待压缩图片
% output: 
% datacodebits：压缩后图片数据

%% 量化表
Y_Table = [16, 11, 10, 16, 24, 40, 51, 61 ;
           12, 12, 14, 19, 26, 58, 60, 55 ;
           14, 13, 16, 24, 40, 57, 69, 56 ;
           14, 17, 22, 29, 51, 87, 80, 62 ;
           18, 22, 37, 56, 68, 109,103,77 ;
           24, 35, 55, 64, 81, 104,113,92 ;
           49, 64, 78, 87, 103,121,120,101;
           72, 92, 95, 98, 112,100,103,99 ];  % 亮度量化表

CbCr_Table = [17, 18, 24, 47, 99, 99, 99, 99 ;
              18, 21, 26, 66, 99, 99, 99, 99 ;
              24, 26, 56, 99, 99, 99, 99, 99 ;
              47, 66, 99 ,99, 99, 99, 99, 99 ;
              99, 99, 99, 99, 99, 99, 99, 99 ;
              99, 99, 99, 99, 99, 99, 99, 99 ;
              99, 99, 99, 99, 99, 99, 99, 99 ;
              99, 99, 99, 99, 99, 99, 99, 99 ];   % 色差量化表
          
%% rgb转换为ycbcr
img_ycbcr = rgb2ycbcr(img); 
[row,col,~] = size(img_ycbcr);

%% 进行图片扩展
row_expand = ceil(row/16)*16; %行数上取整再乘16，及扩展成16的倍数
if mod(row,16) ~= 0           %行数不是16的倍数，用最后一行进行扩展
   for i = row:row_expand
       img_ycbcr(i,:,:) = img_ycbcr(row,:,:);
   end
end

col_expand = ceil(col/16)*16; %列数上取整
if mod(col,16) ~= 0           %列数不是16的倍数，用最后一列进行扩展
   for j = col:col_expand
       img_ycbcr(:,j,:) = img_ycbcr(:,col,:);
   end
end

%% 对Y,Cb,Cr分量进行4:2:0采样
Y = img_ycbcr(:,:,1); % Y分量
Cb = zeros(row_expand/2,col_expand/2); % Cb分量
Cr = zeros(row_expand/2,col_expand/2); % Cr分量

for i = 1:row_expand/2
    for j = 1:2:col_expand/2-1         % 奇数列
       Cb(i,j) = double(img_ycbcr(i*2-1,j*2-1,2));     
       Cr(i,j) = double(img_ycbcr(i*2-1,j*2+1,3));     
    end
end
for i = 1:row_expand/2
    for j = 2:2:col_expand/2           %偶数列
       Cb(i,j) = double(img_ycbcr(i*2-1,j*2-2,2));     
       Cr(i,j) = double(img_ycbcr(i*2-1,j*2,3));     
    end
end

%% 对三个通道进行量化
QuanityFactor = 0.5; % 量化因子
Y_dct_q = Dct_Quantize(Y,QuanityFactor,Y_Table);
Cb_dct_q = Dct_Quantize(Cb,QuanityFactor,CbCr_Table);
Cr_dct_q = Dct_Quantize(Cr,QuanityFactor,CbCr_Table);

data = [Y_dct_q(:);Cb_dct_q(:);Cr_dct_q(:)];
data = data+128;

%% 对数据进行RLC编码
datacode = RLCEncode(data);
% datacode = data;

end