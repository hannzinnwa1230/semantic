function img = imgintegration(img1,img2)
img1 = imread("UESTC.bmp");
img2 = imread("lena.tiff");

height = size(img1,1);
width = size(img1,2)+size(img2,2);
for i = size(img2,1)
    for j = size(img2,2)
        if(img2(i,j) == 0)
            img2new(i,j) = 0;
        else
            img2new(i,j) = 255; 
        end
    end
end
imgtemp = zeros(height,width,3);
imgtemp(1:height,1:size(img1,2),:) = img1;
imgtemp(height-size(img2,1):height,size(img1,2)+1:end,1) = img2;
imgtemp(height-size(img2,1):height,size(img1,2)+1:end,2) = img2;
imgtemp(height-size(img2,1):height,size(img1,2)+1:end,3) = img2;
img = imgtemp;

end