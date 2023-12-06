function PSNR = PSNRCalc(image1,image2)
%% 
% 功能：计算PSNR
% input：
% image1：第一幅图像
% image2：第二幅图像
% output：
% PSNR：计算得到的峰值信噪比

%%
imgsize = size(image1);
PSNR = sum(sum(sum((image1-image2).^2)));
PSNR = 10*log10(255^2*imgsize(1)*imgsize(2)/PSNR);

end