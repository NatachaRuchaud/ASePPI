function dst = pix(img, voisinage1, voisinage2)
size_img = size(img);
 S = size(size_img);
% if(S(2)==3)
%     img = rgb2gray(img);
% end
dst = img;
moyenne1 = 0;
moyenne2 = 0;
moyenne3 = 0;
tmp1 = voisinage1;
tmp2 = voisinage2;
for i=1:voisinage1:size_img(1)-voisinage1
    for j = 1:voisinage2:size_img(2)-voisinage2
%         if(i<=size_img(1)-voisinage1)
%            if(j<=size_img(2)-voisinage2)
        while(i+voisinage1-1>size_img(1))
            voisinage1 = round((voisinage1/2)/4)*4;
        end
        while(j+voisinage2-1>size_img(2))
            voisinage2 = round((voisinage2/2)/4)*4;
        end
                for vi=0:1:voisinage1-1
                    for vj= 0:1:voisinage2-1
                        pixelVal1 = img(i + vi, j + vj, 1);
                        moyenne1 = double(moyenne1) + double(pixelVal1);
                        pixelVal2 = img(i + vi, j + vj, 2);
                        moyenne2 = double(moyenne2) + double(pixelVal2);
                        pixelVal3 = img(i + vi, j + vj, 3);
                        moyenne3 = double(moyenne3) + double(pixelVal3);
                    end
                end
                moyenne1 = round(moyenne1/((voisinage1)*(voisinage2)));
                moyenne2 = round(moyenne2/((voisinage1)*(voisinage2)));
                moyenne3 = round(moyenne3/((voisinage1)*(voisinage2)));
                for vi=0:1:voisinage1-1
                    for vj= 0:1:voisinage2-1
                        %dst(i+vi, j+vj)= moyenne1;
                        dst(i+vi, j+vj, 1)= moyenne1;
                        dst(i+vi, j+vj, 2)= moyenne2;
                        dst(i+vi, j+vj, 3)= moyenne3;
                    end
                end
                moyenne1 = 0;
                moyenne2 = 0;
                moyenne3 = 0;
                voisinage1 = tmp1;
                voisinage2 = tmp2;
            %end
        %end
    end
end
 %rect = [1, 1, size_img(2)-voisinage2, size_img(1)-voisinage1];
 %dst = imcrop(dst, rect);
 %dst = imresize(dst, [size_img(1), size_img(2)]);


%%COMPUTE DIFF
% imwrite(dst, 'protected.jpg', 'jpg');
% 
% clear dst
% dst = imread('protected.jpg');
% diff = double(img)-double(dst);
% minim = min(min(min(diff)));
% maxim = max(max(max(diff)));
% diff = (diff-minim)/(maxim-minim)*255;
% quant = round(diff/10);
% rng(1234)
% quant = uint8(quant)+floor(rand(1)*213);
% imwrite(quant,'diff.jp2','jp2');
% clear quant dst
% dst = imread('protected.jpg');
% quant= imread('diff.jp2');
% 
% rng(1234)
% quant = quant-floor(rand(1)*213);
% iquant = quant*10;
% 
% idiff = double(iquant)/255*(maxim-minim)+minim;
% rec = idiff+double(dst);
% imshow(uint8(rec))

%imwrite(dst, strcat(pathOut, name));
end
