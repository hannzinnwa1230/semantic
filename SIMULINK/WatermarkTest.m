clc;
clear;
close all;

img = imread("lena.tiff");
% img = imresize(img,[400,300]);

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
test=repelem(test,10,10);
% test=255*test./max(max(test));
% test=floor(test);
test_gray=mat2gray(test,[0 255])
watermark=gray2rgb(test_gray);

% watermark = imread("uestc.bmp");
% alpha = AlphaGet(img,watermark);
% alpha1 = 40;
% alpha2 = 20;
% [img_w,realwatermark] = AddWatermark(img,watermark,alpha1,alpha2);

alpha = 0.9;
[img_w,realwatermark] = AddWatermark_PN(img,watermark,alpha);

figure(1)
subplot(221)
imshow(img)
title("原始图像")
subplot(222)
imshow(img_w)
title("添加水印后图像")
PSNR = PSNRCalc(img,img_w)

% 
% txData = DataForm(img_w);
% tximage = img_w;
% img_rec = DataDisForm(txData,tximage);

img_rec = img_w;

watermark_pick = PickWatermark_PN(img_rec,alpha);

% figure(2)
subplot(223)
imshow(realwatermark)
title("水印图像")
subplot(224)
imshow(watermark_pick)
title("提取的水印图像")

NC = NCCalc(realwatermark,watermark_pick)

% txData = DataForm(img);
% tximage = img;
% img_rec1 = DataDisForm(txData,tximage);
% PSNR1 = PSNRCalc(img,img_rec1)
% PSNR2 = PSNRCalc(img,img_rec)