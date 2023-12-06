function Matrix_out = Dct_Quantize(Matrix_in,Qua_Factor,Qua_Table)
%% 
% 功能：Dct量化
% input：
% Matrix_in：输入矩阵
% Qua_Factor: 量化因子
% Qua_Table：量化表
% output：
% Matrix_out：量化后的矩阵
%%
Matrix_in = double(Matrix_in)-128;  %层次移动128个灰度级
Matrix_in = blkproc(Matrix_in,[8 8],'dct2(x)');
Qua_Matrix = Qua_Factor.*Qua_Table;             %量化矩阵
Matrix_in = blkproc(Matrix_in,[8 8],'round(x./P1)',Qua_Matrix); %量化，四舍五入
Matrix_out = Matrix_in;         %得到量化后的矩阵
end