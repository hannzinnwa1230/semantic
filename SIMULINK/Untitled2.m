clc;
clear;
close all;

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
% test=255*test./max(max(test));
% test=floor(test);
test_gray=mat2gray(test,[0 255])
originalwatermark=gray2rgb(test_gray);
originalimage = imread("./TrainPic/Test/lena.tiff");
alpha = 60;%7.7433(44)7.9342(46)8.0478(48)
% originalimage = imread("peppers.tiff");
% alpha = 60.0245;%8.0255(44)8.0085(46)8.0245(48)
% originalimage = imread("baboon.tiff");
% alpha = 61.8893;%7.4966(44)7.6078(46)7.8893(48)
% originalimage = imread("airplane.tiff");
% alpha = 60.9863;%7.9969(44)7.9873(46)7.9863(48)

[img_watermark,realwatermark] = AddWatermark(originalimage,originalwatermark,alpha);
PSNR = PSNRCalc(originalimage,img_watermark)

watermark_pick1 = PickWatermark(img_watermark,alpha);
% NC = NCCalc(realwatermark,watermark_pick);

%% 亮度调节
img_watermark2 = imadd(img_watermark,+4);
watermark_pick2 = PickWatermark(img_watermark2,alpha);

img_watermark3 = imadd(img_watermark,-4);
watermark_pick3 = PickWatermark(img_watermark3,alpha);

%% 噪声攻击
img_watermark4 = imnoise(img_watermark,'gaussian',0,0.01);
img_watermark5 = imnoise(img_watermark,'gaussian',0,0.02);
img_watermark6 = imnoise(img_watermark,'salt & pepper',0.01);
img_watermark7 = imnoise(img_watermark,'salt & pepper',0.03);

watermark_pick4 = PickWatermark(img_watermark4,alpha);
watermark_pick5 = PickWatermark(img_watermark5,alpha);
watermark_pick6 = PickWatermark(img_watermark6,alpha);
watermark_pick7 = PickWatermark(img_watermark7,alpha);

%% 滤波攻击
R=img_watermark(:,:,1);
G=img_watermark(:,:,2);
B=img_watermark(:,:,3);
 
R1=medfilt2(R,[3,3]);
G1=medfilt2(G,[3,3]);
B1=medfilt2(B,[3,3]);

R2=medfilt2(R,[5,5]);
G2=medfilt2(G,[5,5]);
B2=medfilt2(B,[5,5]);

%gaussian average
R3=filter2(fspecial('average',3),R);
G3=filter2(fspecial('average',3),G);
B3=filter2(fspecial('average',3),B);

R4=filter2(fspecial('average',5),R);
G4=filter2(fspecial('average',5),G);
B4=filter2(fspecial('average',5),B);

R5=filter2(fspecial('gaussian',3,0.01),R);
G5=filter2(fspecial('gaussian',3,0.01),G);
B5=filter2(fspecial('gaussian',3,0.01),B);

R6=filter2(fspecial('gaussian',5,0.01),R);
G6=filter2(fspecial('gaussian',5,0.01),G);
B6=filter2(fspecial('gaussian',5,0.01),B);
 
img_watermark8 = cat(3,R1,G1,B1); 
img_watermark9 = cat(3,R2,G2,B2); 
img_watermark10 = cat(3,R3,G3,B3); 
img_watermark11 = cat(3,R4,G4,B4); 
img_watermark12 = cat(3,R5,G5,B5); 
img_watermark13 = cat(3,R5,G5,B5); 

watermark_pick8 = PickWatermark(img_watermark8,alpha);
watermark_pick9 = PickWatermark(img_watermark9,alpha);
watermark_pick10 = PickWatermark(img_watermark10,alpha);
watermark_pick11 = PickWatermark(img_watermark11,alpha);
watermark_pick12 = PickWatermark(img_watermark12,alpha);
watermark_pick13 = PickWatermark(img_watermark13,alpha);

%% 裁剪攻击
img_watermark14 = img_watermark;
img_watermark15 = img_watermark;
img_watermark14(1:128,1:128,:) = 0;
img_watermark15(1:256,1:128,:) = 0;
watermark_pick14 = PickWatermark(img_watermark14,alpha);
watermark_pick15 = PickWatermark(img_watermark15,alpha);

%% 压缩攻击
txData = DataForm(img_watermark);
tximage = img_watermark;
img_rec = DataDisForm(txData,tximage);
img_watermark16 = img_rec;
watermark_pick16 = PickWatermark(img_watermark16,alpha);

% subplot(441)
% imshow(watermark_pick1)
% title('无攻击');
% subplot(442)
% imshow(watermark_pick2)
% title('亮度调节(+4)');
% subplot(443)
% imshow(watermark_pick3)
% title('亮度调节(-4)');
% subplot(444)
% imshow(watermark_pick4)
% title('高斯噪声(σ=0.01)');
% subplot(445)
% imshow(watermark_pick5)
% title('高斯噪声(σ=0.02)');
% subplot(446)
% imshow(watermark_pick6)
% title('均值噪声(σ=0.01)');
% subplot(447)
% imshow(watermark_pick7)
% title('均值噪声(σ=0.03)');
% subplot(448)
% imshow(watermark_pick8)
% title('中值滤波(3×3)');
% subplot(449)
% imshow(watermark_pick9)
% title('中值滤波(5×5)');
% subplot(4,4,10)
% imshow(watermark_pick10)
% title('均值滤波(3×3)');
% subplot(4,4,11)
% imshow(watermark_pick11)
% title('均值滤波(5×5)');
% subplot(4,4,12)
% imshow(watermark_pick12)
% title('高斯滤波(3×3)');
% subplot(4,4,13)
% imshow(watermark_pick13)
% title('高斯滤波(5×5)');
% subplot(4,4,14)
% imshow(watermark_pick14)
% title('裁剪攻击(1/16)');
% subplot(4,4,15)
% imshow(watermark_pick15)
% title('裁剪攻击(1/8)');
% subplot(4,4,16)
% imshow(watermark_pick16)
% title('JPEG压缩');

% subplot(441)
% imshow(watermark_pick1)
% title('无攻击');
subplot(531)
imshow(watermark_pick2)
title('亮度调节(+4)');
subplot(532)
imshow(watermark_pick3)
title('亮度调节(-4)');
subplot(533)
imshow(watermark_pick4)
title('高斯噪声(σ=0.01)');
subplot(534)
imshow(watermark_pick5)
title('高斯噪声(σ=0.02)');
subplot(535)
imshow(watermark_pick6)
title('均值噪声(σ=0.01)');
subplot(536)
imshow(watermark_pick7)
title('均值噪声(σ=0.03)');
subplot(537)
imshow(watermark_pick8)
title('中值滤波(3×3)');
subplot(538)
imshow(watermark_pick9)
title('中值滤波(5×5)');
subplot(539)
imshow(watermark_pick10)
title('均值滤波(3×3)');
subplot(5,3,10)
imshow(watermark_pick11)
title('均值滤波(5×5)');
subplot(5,3,11)
imshow(watermark_pick12)
title('高斯滤波(3×3)');
subplot(5,3,12)
imshow(watermark_pick13)
title('高斯滤波(5×5)');
subplot(5,3,13)
imshow(watermark_pick14)
title('裁剪攻击(1/16)');
subplot(5,3,14)
imshow(watermark_pick15)
title('裁剪攻击(1/8)');
subplot(5,3,15)
imshow(watermark_pick16)
title('JPEG压缩');


NC = NCCalc(realwatermark,watermark_pick1)
NC = NCCalc(realwatermark,watermark_pick2)
NC = NCCalc(realwatermark,watermark_pick3)
NC = NCCalc(realwatermark,watermark_pick4)
NC = NCCalc(realwatermark,watermark_pick5)
NC = NCCalc(realwatermark,watermark_pick6)
NC = NCCalc(realwatermark,watermark_pick7)
NC = NCCalc(realwatermark,watermark_pick8)
NC = NCCalc(realwatermark,watermark_pick9)
NC = NCCalc(realwatermark,watermark_pick10)
NC = NCCalc(realwatermark,watermark_pick11)
NC = NCCalc(realwatermark,watermark_pick12)
NC = NCCalc(realwatermark,watermark_pick13)
NC = NCCalc(realwatermark,watermark_pick14)
NC = NCCalc(realwatermark,watermark_pick15)
NC = NCCalc(realwatermark,watermark_pick16)