clear all; close all; dbstop if error;

% Processing Parameters
methodStr = {'nearest','bilinear','bicubic'};
waveletStr = {'db1'};
iterations = {1, 2, 3, 4};
n_iterations = length(methodStr)*length(waveletStr)*length(iterations);
alpha = 2;

% Input Data
load woman;
img_original = double(X);

% Data Holders
rpsnr_low = zeros(1,n_iterations);
rpsnr_high = zeros(1,n_iterations);
rssim_low = zeros(1,n_iterations);
rssim_high = zeros(1,n_iterations);
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
                [LL,LH,HL,HH] = dwt2(img_original,waveletStr{ii});
                
            case 2
                [LL,~,~,~] = dwt2(img_original,waveletStr{ii});
                [LL,LH,HL,HH] = dwt2(LL,waveletStr{ii});
                
            case 3
                [LL,~,~,~] = dwt2(img_original,waveletStr{ii});
                [LL,~,~,~] = dwt2(LL,waveletStr{ii});
                [LL,LH,HL,HH] = dwt2(LL,waveletStr{ii});
                
            case 4
                [LL,~,~,~] = dwt2(img_original,waveletStr{ii});
                % Add SR inbetween 
                [LL,~,~,~] = dwt2(LL,waveletStr{ii});
                [LL,~,~,~] = dwt2(LL,waveletStr{ii});
                [LL,LH,HL,HH] = dwt2(LL,waveletStr{ii});
                
            otherwise
                error('Unsupported iteration level number');
        end

        for jj=1:length(methodStr)

            idx = idx+1;

            % Perform Interpolation
            jHL = imresize(HL, alpha^iterations{kk}, methodStr{jj});
            jLH = imresize(LH, alpha^iterations{kk}, methodStr{jj});
            jHH = imresize(HH, alpha^iterations{kk}, methodStr{jj});
            jin = imresize(img_original, alpha/2, methodStr{jj});

            % Perform Reconstruction
            %img_rec = idwt2(img_original,jLH,jHL,jHH,waveletStr{ii},size(img_original).*alpha);
            img_rec = idwt2(jin,jLH,jHL,jHH,waveletStr{ii},size(img_original).*alpha);

            % Compute Quantitative Evaluation Parameters
            rec_resized = imresize(img_rec, size(img_original), 'bicubic');
            err = img_original - rec_resized;
            err = err.*err;
            imse = sum(err(:))/(size(img_original,1)*size(img_original,2));
            rpsnr_low(idx) = 10*log10(1*1/imse); % R = 1 for double precision floating-point data type
            rssim_low(idx) = ssim(rec_resized, img_original);

            og_resized = imresize(img_original, size(img_rec), 'bicubic');
            err = og_resized - img_rec;
            err = err.*err;
            imse = sum(err(:))/(size(img_original,1)*size(img_original,2));
            rpsnr_high(idx) = 10*log10(1*1/imse); % R = 1 for double precision floating-point data type
            rssim_high(idx) = ssim(img_rec,og_resized);

            rniqe(idx) = niqe(img_rec);

            % Display Results
            figure; set(gcf, 'Position', dispSize); subplot(2,2,1); imagesc(img_original); colormap(gray);
            title(sprintf('%d x %d Input Image',size(img_original,1),size(img_original,2)));

            subplot(2,2,2); imagesc([LL LH; HL HH]); colormap(gray);
            title(sprintf('%s Level %d Sub-band Decomposition',waveletStr{ii},iterations{kk}));

            subplot(2,2,3:4); imagesc(img_rec); colormap(gray); xlim([-88 600]);
            title(sprintf('%d x %d Level %d Reconstructed Image\n%s Wavelet, %s Interpolation',size(img_rec,1),size(img_rec,2),iterations{kk},waveletStr{ii},methodStr{jj}));
        end
    end
end

figure; set(gcf, 'Position', dispSize); 
subplot(3,1,1); plot((1:1:n_iterations), rpsnr_low); grid on; hold on;
plot((1:1:n_iterations),rpsnr_high); ylabel('PSNR'); xlabel('Case Number');
title('Super-Resolved Images PSNR');
legend('256x256 MSE','512x512 MSE','location','best');

subplot(3,1,2); plot((1:1:n_iterations),rssim_low); grid on; hold on;
plot((1:1:n_iterations),rssim_high); ylabel('SSIM'); xlabel('Case Number');
title('Super-Resolved Images Structural Similarity Index Measure');
legend('256x256','512x512','location','best');

subplot(3,1,3); plot((1:1:n_iterations),rniqe); grid on;
ylabel('NIQE'); xlabel('Case Number');
title(sprintf('Super-Resolved Images Naturalness Image Quality Evaluator\nInput Image Score = %f',niqe(img_original)));