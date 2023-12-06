function img_arnold = arnold(img,num,a,b)
%% 
% 功能：对输入图像进行arnold变换
% input:
% img: 输入图像
% num：变换次数
% a，b：变换系数
% output:
% img_arnold：输出图像

%% 进行arnold变换
IMG = img;
size_img = size(IMG);
h = size_img(1);
for n = 1:num
    for y = 1:h
        for x = 1:h
            %防止取余过程中出现错误，先把坐标系变换成从0 到 N-1
            xx = mod((x-1)+(y-1)*b,h)+1;
            yy = mod((x-1)*a+(y-1)*(a*b+1),h)+1;
            IMG_arnold(yy,xx) = IMG(y,x);
        end
    end
end
img_arnold = IMG_arnold;

end