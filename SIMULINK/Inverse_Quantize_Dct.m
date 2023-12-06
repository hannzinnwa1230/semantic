function Matrix_out = Inverse_Quantize_Dct(Matrix_in,Qua_Factor,Qua_Table)
%% 
% 功能：iDct反量化
% input：
% Matrix_in：输入矩阵
% Qua_Factor: 量化因子
% Qua_Table：量化表
% output：
% Matrix_out：反量化后的矩阵
%% 
Qua_Matrix = Qua_Factor.*Qua_Table;    
Matrix_in = blkproc(Matrix_in,[8 8],'x.*P1',Qua_Matrix);
[row,column] = size(Matrix_in);
Matrix_in = blkproc(Matrix_in,[8 8],'idct2(x)'); 
Matrix_in = uint8(Matrix_in+128);
for i = 1:row
   for j = 1:column
       if Matrix_in(i,j) > 255
           Matrix_in(i,j) = 255;
       elseif Matrix_in(i,j) < 0
           Matrix_in(i,j) = 0;
       end
   end
end
Matrix_out = Matrix_in;      
end