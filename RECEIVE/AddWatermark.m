function [img_watermark,realwatermark] = AddWatermark(originalimage,originalwatermark,alpha)
%% 
% 功能：对图像进行水印添加
% input：
% originalimage：原始图像
% originalwatermark：水印图像
% alpha：水印修正因子
% output：
% img_watermark：添加水印的图像
% realwatermark：实际嵌入的水印图像

%% 对原始图像进行处理
size_oi = size(originalimage);
% height1 = ceil((512-size_oi(1))/2);
% height2 = floor((512-size_oi(1))/2);
% width1 = ceil((512-size_oi(2))/2);
% width2 = floor((512-size_oi(2))/2);
% image = zeros(512,512,3);
% image(height1+1:height1+size_oi(1),width1+1:width1+size_oi(2),:) = originalimage;
image = imresize(originalimage,[512,512]);
YCbCr = rgb2ycbcr(image);

%% 对水印图像进行处理
realwatermark = imresize(originalwatermark,[64,64]);
realwatermark = imbinarize(rgb2gray(realwatermark));
watermark = arnold(realwatermark,10,1,1);
size_wm = size(watermark);

%% 生成A矩阵
img = YCbCr;
[cA1,cH1,cV1,cD1] = dwt2(img,'haar');
img_g = cA1(:,:,1);
img_g_dct = blkproc(img_g,[4,4],'dct2');
B = mat2cell(img_g_dct,4*ones(1,64),4*ones(1,64));
% for i = 1:size_wm(1)/2
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         beta = alpha(1);
%         Z = mod(S(1,1),beta);
%         if(watermark(i,j) == 0)
%             n = round(Z/beta);
%             m = Z/beta;
%             if(n<m)
%                 S(1,1) = S(1,1)-Z+(n+1/4)*beta;
%             else
%                 S(1,1) = S(1,1)-Z+(n-3/4)*beta;
%             end
%         elseif(watermark(i,j) == 1)
%             n = round(Z/beta);
%             m = Z/beta;
%             if(n<m)
%                 S(1,1) = S(1,1)-Z+(n+3/4)*beta;
%             else
%                 S(1,1) = S(1,1)-Z+(n-1/4)*beta;
%             end
%         end        
%         A_w = U*S*V';
%         A_w = mat2cell(A_w,4,4);
%         B_w(i,j) = A_w;
%     end
% end
% for i = size_wm(1)/2+1:size_wm(1)
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         omiga = alpha(2);
%         beta = round(S(1,1)/omiga);
%         if(watermark(i,j) == 0)
%             if(mod(beta,2) == 1)
%                 S(1,1) = (beta+1)*omiga;
%             else
%                 S(1,1) = beta*omiga;
%             end
%         elseif(watermark(i,j) == 1)
%             if(mod(beta,2) == 1)
%                 S(1,1) = beta*omiga;
%             else
%                 S(1,1) = (beta+1)*omiga;
%             end
%         end        
%         A_w = U*S*V';
%         A_w = mat2cell(A_w,4,4);
%         B_w(i,j) = A_w;
%     end
% end



% for i = 1:size_wm(1)
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         beta = alpha*(1/3*(S(2,2)+S(3,3)+S(4,4)));
%         Z = mod(S(1,1),beta);
%         if(watermark(i,j) == 0)
%             S(1,1) = S(1,1)-Z+1/4*beta;
%         elseif(watermark(i,j) == 1)
%             S(1,1) = S(1,1)-Z+3/4*beta;
%         end        
%         A_w = U*S*V';
%         A_w = mat2cell(A_w,4,4);
%         B_w(i,j) = A_w;
%     end
% end
% for i = 1:size_wm(1)
%     for j = 1:size_wm(2)
%         A = cell2mat(B(i,j));
%         [U,S,V] = svd(A);
%         if(watermark(i,j) == 0)
%             S(2,2) = (1-alpha)*S(2,2)+alpha*S(3,3);
%         elseif(watermark(i,j) == 1)
%             S(2,2) = (1-alpha)*S(2,2)+alpha*S(4,4);
%         end        
%         A_w = U*S*V';
%         A_w = mat2cell(A_w,4,4);
%         B_w(i,j) = A_w;
%     end
% end

for i = 1:size_wm(1)
    for j = 1:size_wm(2)
        A = cell2mat(B(i,j));
        [U,S,V] = svd(A);
        omiga = alpha;
        beta = round(S(1,1)/omiga);
        if(watermark(i,j) == 0)
            if(mod(beta,2) == 1)
                S(1,1) = (beta+1)*omiga;
            else
                S(1,1) = beta*omiga;
            end
        elseif(watermark(i,j) == 1)
            if(mod(beta,2) == 1)
                S(1,1) = beta*omiga;
            else
                S(1,1) = (beta+1)*omiga;
            end
        end        
        A_w = U*S*V';
        A_w = mat2cell(A_w,4,4);
        B_w(i,j) = A_w;
    end
end




B_w = cell2mat(B_w);

%% 重建图像
image_gw_dct = B_w;
image_gw = blkproc(image_gw_dct,[4,4],'idct2');
cA1(:,:,1) = image_gw;
image_w = idwt2(cA1,cH1,cV1,cD1,'haar');
image_w = ycbcr2rgb(uint8(image_w));

%% 得到添加水印的图像
% img_watermark = image_w(height1+1:height1+size_oi(1),width1+1:width1+size_oi(2),:);
img_watermark = (imresize(image_w,[size_oi(1),size_oi(2)]));

end
