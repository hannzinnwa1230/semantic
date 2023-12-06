function watermark_pick = PickWatermark(image,alpha)
%% 
% 功能：对图像进行水印提取
% input：
% image：原始图像
% alpha：水印修正因子
% U_key,S_key,V_key：密匙矩阵
% output：
% iwatermark：提取的水印图像

%% 处理原始图像
image_rec = imresize(image,[512,512]);
size_img = size(image_rec);
size_wm = [64,64];
YCbCr = rgb2ycbcr(image_rec);

%% 提取水印
img = YCbCr;
[cA1,cH1,cV1,cD1] = dwt2(img,'haar');
image_pick = cA1(:,:,1);
image_pick_dct = blkproc(image_pick,[4,4],'dct2');
B = mat2cell(image_pick_dct,4*ones(1,64),4*ones(1,64));
% for i = 1:size_wm(1)/2
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         beta = alpha(1);
%         Z = mod(S(1,1),beta);
%         if(Z <= beta/2)
%             watermark_pick(i,j) = 0;
%         else
%             watermark_pick(i,j) = 1;
%         end
%     end
% end
% for i = size_wm(1)/2+1:size_wm(1)
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         omiga = alpha(2);
%         beta = round(S(1,1)/omiga);
%         if(mod(beta,2) == 0)
%             watermark_pick(i,j) = 0;
%         else
%             watermark_pick(i,j) = 1;
%         end
%     end
% end

for i = 1:size_wm(1)
    for j = 1:size_wm(2)
        A = cell2mat(B(i,j));
        [U,S,V] = svd(A);
        omiga = alpha;
        beta = round(S(1,1)/omiga);
        if(mod(beta,2) == 0)
            watermark_pick(i,j) = 0;
        else
            watermark_pick(i,j) = 1;
        end
    end
end

watermark_pick = logical(arnoldrec(watermark_pick,10,1,1));

end
