function NC = NCCalc(image1,image2)
%%
% 功能：计算两幅图片相关性
% input：
% image1：第一幅图片
% image2：第二幅图片
% output：
% NC：计算得到的相关性

%%
NC1 = sum(sum(image1.*image2));
NC2 = sum(sum(image1.*image1));
NC = min(NC1/NC2,NC2/NC1);

end