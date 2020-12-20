% Input Data
%load woman;
X = imread('lena_std.tif');
X = rgb2gray(X);
img_original = double(X);
alpha = 2;

waveStr = 'db1';
intStr = 'bicubic';

% Decomposition
[LL,LH,HL,HH] = dwt2(img_original,waveStr);

% Interpolation
jHL = imresize(HL, alpha, intStr);
jLH = imresize(LH, alpha, intStr);
jHH = imresize(HH, alpha, intStr);
im = imresize(img_original, alpha/2, intStr);

% Reconstruction
img_rec = idwt2(im,jLH,jHL,jHH,waveStr,size(img_original).*alpha).*alpha;

% Display Results
d = imresize(img_original, alpha, intStr);
figure; subplot(1,2,1); imagesc(img_original); colormap(gray);
title(sprintf('%d x %d Input Image',size(img_original,1),size(img_original,2)));

subplot(1,2,2); imagesc(img_rec); colormap(gray); %xlim([-50 size(img_rec,2)+50]);
title(sprintf('%d x %d Level %d Reconstructed Image\n%s Wavelet, %s Interpolation',...
size(img_rec,1),size(img_rec,2),1,waveStr,intStr));

fprintf('upsampled original PSNR: %f\n',psnr(imresize(img_original, alpha, intStr),...
                                             img_rec));
fprintf('upsampled original SSIM: %f\n',ssim(imresize(img_original, alpha, intStr),...
                                             img_rec));
fprintf('downsampled result PSNR: %f\n',psnr(img_original,...
                                             imresize(img_rec, 1/alpha, intStr)));
fprintf('downsampled result SSIM: %f\n',ssim(img_original,...
                                             imresize(img_rec, 1/alpha, intStr)));
