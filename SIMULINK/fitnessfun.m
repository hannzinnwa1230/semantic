function output = fitnessfun(originalimage,originalwatermark,alpha,beta)

size_alpha = size(alpha);
for i = 1:size_alpha(1)
    [img_watermark,realwatermark] = AddWatermark(originalimage,originalwatermark,alpha(i,1));
    PSNR = PSNRCalc(originalimage,img_watermark);
    watermark_pick = PickWatermark(img_watermark,alpha(i,1));
    NC = NCCalc(realwatermark,watermark_pick);
    output(i,1) = PSNR+beta*NC;
end

end