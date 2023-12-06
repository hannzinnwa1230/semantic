clc;
clear;
close all;
% 
% load alpha;
% originalimage = imread("lena.bmp");
% originalwatermark = imread("uestc.bmp");
% size_oi = size(originalimage);
% image = imresize(originalimage,[512,512]);
% YCbCr = rgb2ycbcr(image);

% %% 对水印图像进行处理
% realwatermark = imresize(originalwatermark,[64,64]);
% realwatermark = imbinarize(rgb2gray(realwatermark));
% watermark = arnold(realwatermark,10,1,1);
% size_wm = size(watermark);
% 
% %% 生成A矩阵
% img = YCbCr;
% [cA1,cH1,cV1,cD1] = dwt2(img,'haar');
% img_g = cA1(:,:,1);
% img_g_dct = blkproc(img_g,[4,4],'dct2');
% B = mat2cell(img_g_dct,4*ones(1,64),4*ones(1,64));
% %% 嵌入水印
% for i = 1:size_wm(1)
%     for j = 1:size_wm(2)
%             A = cell2mat(B(i,j));
%             [U,S,V] = svd(A);
%             beta = alpha(i,j);
%             Z = mod(S(1,1),beta);
%             if(watermark(i,j) == 0)
%                 n = round(Z/beta);
%                 m = Z/beta;
%                 if(n<m)
%                     S(1,1) = S(1,1)-Z+(n+1/4)*beta;
%                 else
%                     S(1,1) = S(1,1)-Z+(n-3/4)*beta;
%                 end
%             elseif(watermark(i,j) == 1)
%                 n = round(Z/beta);
%                 m = Z/beta;
%                 if(n<m)
%                     S(1,1) = S(1,1)-Z+(n+3/4)*beta;
%                 else
%                     S(1,1) = S(1,1)-Z+(n-1/4)*beta;
%                 end
%             end
%         A_w = U*S*V';
%         A_w = mat2cell(A_w,4,4);
%         B_w(i,j) = A_w;
%     end
% end
% 
% B_w = cell2mat(B_w);
% image_gw_dct = B_w;
% clear B_w;
% image_gw = blkproc(image_gw_dct,[4,4],'idct2');
% cA1(:,:,1) = image_gw;
% image_w = idwt2(cA1,cH1,cV1,cD1,'haar');
% image_w = ycbcr2rgb(uint8(image_w));
% img_watermark = (imresize(image_w,[size_oi(1),size_oi(2)]));
% 
% %% 提取水印
% image_rec = imresize(img_watermark,[512,512]);
% size_img = size(image_rec);
% size_wm = [64,64];
% YCbCr = rgb2ycbcr(image_rec);
% img = YCbCr;
% [cA1,cH1,cV1,cD1] = dwt2(img,'haar');
% image_pick = cA1(:,:,1);
% image_pick_dct = blkproc(image_pick,[4,4],'dct2');
% B = mat2cell(image_pick_dct,4*ones(1,64),4*ones(1,64));
% for i = 1:size_wm(1)
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         beta = alpha(i,j);
%         Z = mod(S(1,1),beta);
%         if(Z <= beta/2)
%             watermark_pick(i,j) = 0;
%         else
%             watermark_pick(i,j) = 1;
%         end
%     end
% end
% 
% PSNR = PSNRCalc(originalimage,img_watermark)
% NC = NCCalc(realwatermark,watermark_pick)

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
test_gray=mat2gray(test,[0 255]);
originalwatermark=gray2rgb(test_gray);

img = imread("./TrainPic/Test/lena.tiff");
watermark = originalwatermark;
ind = 1;

for x = 50:-1:20
    [img_w,realwatermark] = AddWatermark(img,watermark,x);
    watermark_pick = PickWatermark(img_w,x);
    NC = NCCalc(realwatermark,watermark_pick);
    if(NC == 1)
        alpha = x
    else
        break;
    end
end

% alpha = (1:0.5:50);
% Mul = NC.*PSNR;
% index = find(Mul == max(Mul))