clear all; close all; dbstop if error;

% Processing Parameters
methodStr = {'nearest','bilinear','bicubic'};
%methodStr = {'bilinear'};
waveletStr = {'db1','db2','bior1.3','bior2.2'};%, ''};
iterations = {1, 2, 3, 4};
n_iterations = length(methodStr)*length(waveletStr)*length(iterations);
alpha = 2;

% Input Data
%load woman;
%X = checkerboard(32);
%img_original = double(X);
X = imread('lena_std.tif');
X = rgb2gray(X);
img_original = double(X);

% Data Holders
rpsnr = zeros(1,n_iterations);
rssim = zeros(1,n_iterations);
rniqe = zeros(1,n_iterations);

% Display Parameters
scrnSize = get(0, 'screensize');
dispSize = [1 1 scrnSize(3)/2 scrnSize(4)];
idx = 0;

for ii=1:length(waveletStr)
    
    for kk = 1:length(iterations)
    
        % Perform Decomposition per iteration level
        switch iterations{kk}
            
            case 1
                for jj=1:length(methodStr)
                    idx = idx+1;
                    
                    % --- Level 1 --- %
                    % Decomposition
                    img = imresize(img_original, 1/alpha, methodStr{jj});
                    [LL,LH,HL,HH] = dwt2(img,waveletStr{ii});
                    
                    % Interpolation
                    jHL = imresize(HL, alpha, methodStr{jj});
                    jLH = imresize(LH, alpha, methodStr{jj});
                    jHH = imresize(HH, alpha, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img,jLH,jHL,jHH,waveletStr{ii},size(img_original));
                
                    
                    dispBands = [LL LH; HL HH];
                    [rpsnr(idx), rssim(idx), rniqe(idx)] = getStats(img_original, img_rec, img, ...
                                                            dispBands, dispSize, ...
                                                            waveletStr{ii}, iterations{kk}, methodStr{jj});
                end
                
            case 2
                for jj=1:length(methodStr)
                    idx = idx+1;
                
                    % --- Level 1 --- %
                    % Deconstruction
                    img = imresize(img_original, 1/alpha, methodStr{jj});
                    [~,LH,HL,HH] = dwt2(img,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(HL, alpha, methodStr{jj}); % 128
                    jLH = imresize(LH, alpha, methodStr{jj});
                    jHH = imresize(HH, alpha, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img,jLH,jHL,jHH,waveletStr{ii},size(img_original)); %256
                    
                    
                    % --- Level 2 --- %
                    % Deconstruction
                    img_rec = imresize(img_rec,1/alpha,methodStr{jj}); %128
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(LLHL, alpha^2, methodStr{jj});
                    jLH = imresize(LLLH, alpha^2, methodStr{jj});
                    jHH = imresize(LLHH, alpha^2, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img_rec,jLH,jHL,jHH,waveletStr{ii},size(img_original)); % 256
                    
                    dispBands = [LLLL,LLLH;LLHL,LLHH];
                    [rpsnr(idx), rssim(idx), rniqe(idx)] = getStats(img_original, img_rec, img, ...
                                                            dispBands, dispSize, ...
                                                            waveletStr{ii}, iterations{kk}, methodStr{jj});
                end
                
            case 3
                for jj=1:length(methodStr)
                    idx = idx+1;
                
                    % --- Level 1 --- %
                    % Deconstruction
                    img = imresize(img_original, 1/alpha, methodStr{jj});
                    [~,LH,HL,HH] = dwt2(img,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(HL, alpha, methodStr{jj}); % 128
                    jLH = imresize(LH, alpha, methodStr{jj});
                    jHH = imresize(HH, alpha, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img,jLH,jHL,jHH,waveletStr{ii},size(img_original)); %256
                    
                    
                    % --- Level 2 --- %
                    % Deconstruction
                    img_rec = imresize(img_rec,1/alpha,methodStr{jj}); %128
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(LLHL, alpha^2, methodStr{jj});
                    jLH = imresize(LLLH, alpha^2, methodStr{jj});
                    jHH = imresize(LLHH, alpha^2, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img_rec,jLH,jHL,jHH,waveletStr{ii},size(img_original)); % 256
                    
                    
                    % --- Level 3 --- %
                    % Deconstruction
                    img_rec = imresize(img_rec,1/alpha,methodStr{jj}); %128
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(LLHL, alpha^2, methodStr{jj});
                    jLH = imresize(LLLH, alpha^2, methodStr{jj});
                    jHH = imresize(LLHH, alpha^2, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img_rec,jLH,jHL,jHH,waveletStr{ii},size(img_original)); % 256
                
                    dispBands = [LLLL,LLLH;LLHL,LLHH];
                    [rpsnr(idx), rssim(idx), rniqe(idx)] = getStats(img_original, img_rec, img, ...
                                                            dispBands, dispSize, ...
                                                            waveletStr{ii}, iterations{kk}, methodStr{jj});
                end
                
            case 4
                for jj=1:length(methodStr)
                    idx = idx+1;
                
                    % --- Level 1 --- %
                    % Deconstruction
                    img = imresize(img_original, 1/alpha, methodStr{jj});
                    [~,LH,HL,HH] = dwt2(img,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(HL, alpha, methodStr{jj}); % 128
                    jLH = imresize(LH, alpha, methodStr{jj});
                    jHH = imresize(HH, alpha, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img,jLH,jHL,jHH,waveletStr{ii},size(img_original)); %256
                    
                    
                    % --- Level 2 --- %
                    % Deconstruction
                    img_rec = imresize(img_rec,1/alpha,methodStr{jj}); %128
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(LLHL, alpha^2, methodStr{jj});
                    jLH = imresize(LLLH, alpha^2, methodStr{jj});
                    jHH = imresize(LLHH, alpha^2, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img_rec,jLH,jHL,jHH,waveletStr{ii},size(img_original)); % 256
                    
                    
                    % --- Level 3 --- %
                    % Deconstruction
                    img_rec = imresize(img_rec,1/alpha,methodStr{jj}); %128
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(LLHL, alpha^2, methodStr{jj});
                    jLH = imresize(LLLH, alpha^2, methodStr{jj});
                    jHH = imresize(LLHH, alpha^2, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img_rec,jLH,jHL,jHH,waveletStr{ii},size(img_original)); % 256
                    
                    
                    % --- Level 4 --- %
                    % Deconstruction
                    img_rec = imresize(img_rec,1/alpha,methodStr{jj}); %128
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii}); % 64
                    
                    % Interpolation
                    jHL = imresize(LLHL, alpha^2, methodStr{jj});
                    jLH = imresize(LLLH, alpha^2, methodStr{jj});
                    jHH = imresize(LLHH, alpha^2, methodStr{jj});

                    % Reconstruction
                    img_rec = idwt2(img_rec,jLH,jHL,jHH,waveletStr{ii},size(img_original)); % 256
                
                    dispBands = [LLLL,LLLH;LLHL,LLHH];
                    [rpsnr(idx), rssim(idx), rniqe(idx)] = getStats(img_original, img_rec, img, ...
                                                            dispBands, dispSize, ...
                                                            waveletStr{ii}, iterations{kk}, methodStr{jj});
                end
                
            otherwise
                error('Unsupported iteration level number');
        end
    end
end

figure; set(gcf, 'Position', dispSize); 
subplot(2,1,1); plot((1:1:n_iterations), rpsnr); grid on;
ylabel('PSNR'); xlabel('Case Number');
title('Super-Resolved Images PSNR');

subplot(2,1,2); plot((1:1:n_iterations),rssim); grid on;
ylabel('SSIM'); xlabel('Case Number');
title('Super-Resolved Images Structural Similarity Index Measure');

%subplot(3,1,3); plot((1:1:n_iterations),rniqe); grid on;
%ylabel('NIQE'); xlabel('Case Number');
%title(sprintf('Super-Resolved Images Naturalness Image Quality Evaluator\nInput Image Score = %f',niqe(img)));

function [psnr_out, ssim_out, niqe_out] = getStats(img_original, img_reconstructed, img, dispBands, dispSize, waveletStr, levelStr, intStr)

% Compute Quantitative Evaluation Parameters
err = img_original - img_reconstructed;
err = err.*err;
imse = sum(err(:))/(size(img_original,1)*size(img_original,2));
psnr_out = 10*log10(1*1/imse); % R = 1 for double precision floating-point data type
psnr_out = psnr(img_reconstructed, img_original);
ssim_out = ssim(img_reconstructed, img_original);
niqe_out = niqe(img_reconstructed);

% Display Results
figure; set(gcf, 'Position', dispSize); subplot(2,2,1); imagesc(img); colormap(gray);
title(sprintf('%d x %d Input Image',size(img,1),size(img,2)));

subplot(2,2,2); imagesc(dispBands); colormap(gray);
title(sprintf('%s Level %d Sub-band Decomposition',waveletStr,levelStr));

subplot(2,2,3:4); imagesc(img_reconstructed); colormap(gray); xlim([-94 350]);
title(sprintf('%d x %d Level %d Reconstructed Image\n%s Wavelet, %s Interpolation',...
    size(img_reconstructed,1),size(img_reconstructed,2),levelStr,waveletStr,intStr));
end