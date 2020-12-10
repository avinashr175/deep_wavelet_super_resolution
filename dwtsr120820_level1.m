clear all; close all; dbstop if error;

methodStr = {'bicubic'};
waveletStr = {'db1'};%,'db2','bior1.3','bior2.2'};%, ''};
alpha = 2;

X = imread('lena_std.tif');
img_original = rgb2gray(X);
%img_original = double(X);

% Level 1
[LL,LH,HL,HH] = dwt2(img_original,waveletStr{1});

jHL = imresize(HL, alpha, methodStr{1});
jLH = imresize(LH, alpha, methodStr{1});
jHH = imresize(HH, alpha, methodStr{1});

img_rec1 = idwt2(img_original,jLH,jHL,jHH,waveletStr{1},size(img_original).*alpha);
img_rec1 = uint8(rescale(img_rec1,0,255));

figure; subplot(1,2,1); imagesc(img_original); title('Original Image');
subplot(1,2,2); imagesc(img_rec1); title('Level 1 Reconstruction');

figure; subplot(1,2,1); imagesc(img_original); title('Original Image');
img_recd1 = imresize(img_rec1,1/alpha,'bicubic');
subplot(1,2,2); imagesc(img_recd1); 
p = psnr(img_original, img_recd1);
title(sprintf('Downsampled Level 1 Reconstruction\nPSNR = %d',p));

in = imresize(img_original,alpha,'bicubic');
figure; subplot(1,2,1); imagesc(in); title('Upsampled Original Image');
subplot(1,2,2); imagesc(img_rec1); 
p = psnr(in, img_rec1);
title(sprintf('Level 1 Reconstruction\nPSNR = %d',p));