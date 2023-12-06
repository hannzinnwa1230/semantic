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
test_gray=mat2gray(test,[0 255])
originalwatermark=gray2rgb(test_gray);

% originalwatermark = imread("UESTC.bmp");
originalimage = imread("lena.tiff");
alpha = 25;
[img_watermark,realwatermark] = AddWatermark(originalimage,originalwatermark,alpha);
PSNR = PSNRCalc(originalimage,img_watermark)
txData = DataForm(img_watermark);
tximage = img_watermark;
img_rec = DataDisForm(txData,tximage);
img_watermark = img_rec;

% logo = imresize(imread("pic10.jpg"),[64,64]);
% img_logo = img_watermark;
% img_logo(1:64,1:64,:) = logo(1:64,1:64,:);
% imshow(img_logo)

img_logo = img_watermark;

watermark_pick = PickWatermark(img_logo,alpha);
NC = NCCalc(realwatermark,watermark_pick)
