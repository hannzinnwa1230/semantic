clc;
clear;
close all;

%% 参数设置
originalimage = imread("./TrainPic/Test/lena.tiff");
% peppers.tiff
% lena.tiff
% baboon.tiff
% airplane.tiff
% originalwatermark = imread("uestc.bmp");
test=[[15.394878   7.131727   4.3163834 12.776562 ]
 [ 8.295607   6.9973907  6.122822   4.018806 ]
 [ 9.147877   2.8270926  0.9408922  4.7171144]
 [ 3.3111987 11.998605   7.1071925  9.390033 ]
 [ 6.4450126  6.6779184  9.098176   4.424156 ]
 [ 7.128175   3.3501606  7.008772   4.7983623]
 [ 6.9007344  1.246918   5.398112   8.289003 ]
 [ 7.7740874  4.0068054  8.688012   8.622568 ]
 [ 5.8230977 10.019472   8.00594    7.135023 ]
 [13.164159  10.632762  13.188708   5.451745 ]
 [ 9.977087   0.        11.589674   3.6032958]
 [ 1.1209099  5.0947742 11.903207  19.301378 ]
 [ 3.280991  10.447733   9.522085   6.3052506]
 [11.139192  10.568192   7.16008    8.343462 ]
 [ 3.6183174 12.390406   5.469112   2.0599113]
 [ 8.312079   7.8961887 10.188387   6.716921 ]];
test=repelem(test,10,10)
test_gray=mat2gray(test,[0 255])
originalwatermark=gray2rgb(test_gray);
beta = 48;

    %% 定义遗传算法参数
NIND=25;        %个体数目
MAXGEN=100;      %最大遗传代数
PRECI=16;       %变量的二进制位数
GGAP=0.95;      %代沟
px=0.7;         %交叉概率
pm=0.01;        %变异概率
N = 1;          %元素个数
trace=zeros(N+1,MAXGEN);                        %寻优结果的初始值

FieldD=[repmat(PRECI,1,N);repmat([5;10],1,N);repmat([1;0;1;1],1,N)];                      %区域描述器
Chrom=crtbp(NIND,PRECI*N);                      %初始种群
%% 优化
gen=0;                                 %代计数器
X=bs2rv(Chrom,FieldD);                 %计算初始种群的十进制转换
ObjV=fitnessfun(originalimage,originalwatermark,X,beta);       %计算目标函数值
while gen<MAXGEN
   fprintf('%d\n',gen)
   FitnV=ranking(-ObjV);                              %分配适应度值
   SelCh=select('sus',Chrom,FitnV,GGAP);              %选择
   SelCh=recombin('xovsp',SelCh,px);                  %重组
   SelCh=mut(SelCh,pm);                               %变异
   X=bs2rv(SelCh,FieldD);               %子代个体的十进制转换
   ObjVSel=fitnessfun(originalimage,originalwatermark,X,beta);             %计算子代的目标函数值
   [Chrom,ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel); %重插入子代到父代，得到新种群
   X=bs2rv(Chrom,FieldD);
   gen=gen+1;                                             %代计数器增加
   %获取每代的最优解及其序号，Y为最优解,I为个体的序号
   [Y,I]=max(ObjV);
   trace(1:N,gen)=X(I,:);                       %记下每代的最优值
   trace(end,gen)=Y;                               %记下每代的最优值
end
%% 画进化图
plot(1:MAXGEN,trace(end,:));
grid on
xlabel('遗传代数')
ylabel('性能变化')
title('进化过程')
bestY = trace(end,end);
bestX=trace(1:end-1,end);
fprintf(['最优嵌入因子:\nX=',num2str(bestX'),'\n最优性能err=',num2str(bestY),'\n'])

function output = fitnessfun(originalimage,originalwatermark,alpha,beta)

size_alpha = size(alpha);
for i = 1:size_alpha(1)
    [img_watermark,realwatermark] = AddWatermark(originalimage,originalwatermark,alpha(i,:));
    PSNR = PSNRCalc(originalimage,img_watermark)
    watermark_pick = PickWatermark(img_watermark,alpha(i,:));
    NC = NCCalc(realwatermark,watermark_pick)
    output(i,1) = PSNR+beta*NC;
end

end
