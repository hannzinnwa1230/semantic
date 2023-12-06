function img_decode = JPEGDecode(datacode,row,col)
%%
% 功能：对图片进行JPEG解压
% input：
% img：待解压图片数据
% row：原图片的行数
% col：原图片的列数
% output: 
% img_decode：解压后图片

%% 对数据进行RLC解码
data = RLCDecode(datacode);
data = data-128;

%% 计算图片扩展后的尺寸
row_expand = ceil(row/16)*16; 
col_expand = ceil(col/16)*16;

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
          
%% 对三个通道分别反量化和反DCT
% Y_dct_q = reshape(data(1:512*512),[512,512]);
% Cb_dct_q = reshape(data(512*512+1:512*512+256*256),[256,256]);
% Cr_dct_q = reshape(data(512*512+256*256+1:512*512+256*256+256*256),[256,256]);
Y_dct_q = reshape(data(1:row_expand*col_expand),[row_expand,col_expand]);
Cb_dct_q = reshape(data(row_expand*col_expand+1:row_expand*col_expand+row_expand/2*col_expand/2),...
[row_expand/2,col_expand/2]);
Cr_dct_q = reshape(data(row_expand*col_expand+row_expand/2*col_expand/2+1:row_expand*col_expand+...
row_expand/2*col_expand/2+row_expand/2*col_expand/2),[row_expand/2,col_expand/2]);

QuanityFactor = 0.5; % 量化因子
Y_in_q_dct = Inverse_Quantize_Dct(Y_dct_q,QuanityFactor,Y_Table);
Cb_in_q_dct = Inverse_Quantize_Dct(Cb_dct_q,QuanityFactor,CbCr_Table);
Cr_in_q_dct = Inverse_Quantize_Dct(Cr_dct_q,QuanityFactor,CbCr_Table);

%% 恢复出YCBCR图像
YCbCr(:,:,1)=Y_in_q_dct;
for i=1:row_expand/2
   for j=1:col_expand/2
       YCbCr(2*i-1,2*j-1,2)=Cb_in_q_dct(i,j);
       YCbCr(2*i-1,2*j,2)=Cb_in_q_dct(i,j);
       YCbCr(2*i,2*j-1,2)=Cb_in_q_dct(i,j);
       YCbCr(2*i,2*j,2)=Cb_in_q_dct(i,j);  
       YCbCr(2*i-1,2*j-1,3)=Cr_in_q_dct(i,j);
       YCbCr(2*i-1,2*j,3)=Cr_in_q_dct(i,j);
       YCbCr(2*i,2*j-1,3)=Cr_in_q_dct(i,j);
       YCbCr(2*i,2*j,3)=Cr_in_q_dct(i,j);
   end
end

img_decode = ycbcr2rgb(YCbCr);
img_decode(row+1:row_expand,:,:)=[];%去掉扩展的行
img_decode(:,col+1:col_expand,:)=[];%去掉扩展的列

end