clc;
clear;
close all;

% %%
Vec = [];
Tag = [];
alpha = 0.9;
originalwatermark = imread("UESTC.bmp");
realwatermark = imresize(originalwatermark,[64,64]);
realwatermark = imbinarize(rgb2gray(realwatermark));
watermark = arnold(realwatermark,10,1,1);

%% 提取训练数据
for i = 1:10 
    switch i
        case 1        
        originalimage = imread("./TrainPic/pic1.bmp");
        case 2
        originalimage = imread("./TrainPic/pic2.jpg");
        case 3
        originalimage = imread("./TrainPic/lena.bmp");
        case 4
        originalimage = imread("./TrainPic/DZKD.jpg");
        case 5
        originalimage = imread("./TrainPic/pic3.jpg");
        case 6
        originalimage = imread("./TrainPic/pic4.jpg");
        case 7
        originalimage = imread("./TrainPic/pic5.jpg");
        case 8
        originalimage = imread("./TrainPic/pic6.jpg");
        case 9
        originalimage = imread("./TrainPic/pic7.jpg");
        case 10
        originalimage = imread("./TrainPic/pic8.jpg");
        case 11
        originalimage = imread("./TrainPic/pic9.jpg");
        case 12
        originalimage = imread("./TrainPic/pic10.jpg");
        case 13
        originalimage = imread("./TrainPic/pic11.png");
    end
    [image_watermark,~] = AddWatermark_PN(originalimage,originalwatermark,alpha);
    for j = 1:3
        switch j
            case 1        
            image_rec = imadd(image_watermark,2);
            case 2
            image_rec = imnoise(image_watermark,'gaussian',0,0.02);
            case 3
            image_rec = histeq(image_watermark);
        end
        image_rec = imresize(image_rec,[512,512]);
        size_img = size(image_rec);
        size_wm = [64,64];
        YCbCr = rgb2ycbcr(image_rec);
        img = YCbCr;
        [cA1,cH1,cV1,cD1] = dwt2(img,'haar');
        image_pick = cA1(:,:,1);
        image_pick_dct = blkproc(image_pick,[4,4],'dct2');
        B = mat2cell(image_pick_dct,4*ones(1,64),4*ones(1,64));
        for i = 1:size_wm(1)
            for j = 1:size_wm(2)
                A = cell2mat(B(i,j));
                [U,S,V] = svd(A);
                Vectemp((i-1)*size_wm(2)+j,:) = [S(1,1),S(2,2),S(3,3),S(4,4)];
                Tagtemp((i-1)*size_wm(2)+j,:) = (watermark(i,j)+1);
            end
        end     
        Vec = [Vec;Vectemp];
        Tag = [Tag;Tagtemp];
    end
end
save Vec
save Tag


% %% 添加水印
% originalimage = imread("test.jpg");
% img_watermark = AddWatermark(originalimage,originalwatermark,alpha);
% imshow(img_watermark)
% 
% PSNR = PSNRCalc(originalimage,img_watermark)
% 
% image_rec = imresize(img_watermark,[512,512]);
% size_img = size(image_rec);
% size_wm = [64,64];
% YCbCr = rgb2ycbcr(image_rec);
% 
% %% 提取水印
% img = YCbCr;
% [cA1,cH1,cV1,cD1] = dwt2(img,'haar');
% image_pick = cA1(:,:,1);
% image_pick_dct = blkproc(image_pick,[4,4],'dct2');
% B = mat2cell(image_pick_dct,4*ones(1,64),4*ones(1,64));
% Vectest = [];
% for i = 1:size_wm(1)
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         Vectest((i-1)*size_wm(2)+j,:) = [S(1,1),S(2,2),S(3,3),S(4,4),alpha];
%     end
% end
% % load net
% load trainedModel
% % Tagtest = round(sim(net,Vectest')')-1;
% Tagtest = trainedModel.predictFcn(Vectest)-1;
% watermark_pick = reshape((Tagtest),[size_wm(1),size_wm(2)])';
% watermark_pick = logical(arnoldrec(watermark_pick,10,1,1));
% imshow(watermark_pick)
% NC = NCCalc(realwatermark,watermark_pick)


