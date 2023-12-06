function  alpha = AlphaGet(originalimage,originalwatermark)
%% 
% 功能：求得图像最佳嵌入因子
% input：
% originalimage：原始图像
% originalwatermark：水印图像
% output：
% alpha：最佳水印修正因子

NCmax = 0;
alpha = 60;
for x = 60:-1:20
    [img_w,realwatermark] = AddWatermark(originalimage,originalwatermark,x);
    watermark_pick = PickWatermark(img_w,x);
    NC = NCCalc(realwatermark,watermark_pick);
    if(NC >= NCmax)
        NCmax = NC;
        alphatemp = x;
        if(alphatemp < alpha)
            alpha = alphatemp;
        end
    else
        break;
    end
end


end