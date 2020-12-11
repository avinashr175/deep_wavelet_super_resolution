clear all; close all; dbstop if error;

% Vignette Variations
waveletStr = {'db1',... % Haar/db1; asymmetric, orthogonal, biorthogonal, linear phase
              'db2',... % db2; asymmetric, orthogonal, biorthogonal
              'db3',... % db3; asymmetric, orthogonal, biorthogonal
              'db4',... % db4; asymmetric, orthogonal, biorthogonal
              'sym2',... % sym2; near symmetric, orthogonal, biorthogonal (least asymmetric ~ nearly linear phase)
              'sym3',... % sym3; near symmetric, orthogonal, biorthogonal (least asymmetric ~ nearly linear phase)
              'sym4',... % sym4; near symmetric, orthogonal, biorthogonal (least asymmetric ~ nearly linear phase)
              'coif1', ... % coif1; near symmetric, orthogonal, biorthogonal (least asymmetric ~ nearly linear phase)
              'coif2',... % coif2; near symmetric, orthogonal, biorthogonal (least asymmetric ~ nearly linear phase)
              'coif3',... % coif3; near symmetric, orthogonal, biorthogonal (least asymmetric ~ nearly linear phase)
              'bior1.3',... % bior1.3; symmetric, not orthogonal, biorthogonal, linear phase
              'bior1.5',... % bior1.3; symmetric, not orthogonal, biorthogonal, linear phase
              'bior2.2',... % bior2.2; symmetric, not orthogonal, biorthogonal, linear phase
              'bior3.3',... % bior3.3; symmetric, not orthogonal, biorthogonal, linear phase
              'bior3.5',... % bior3.5; symmetric, not orthogonal, biorthogonal, linear phase
              'dmey'}; % discrete meyer; symmetric, not orthogonal, biorthogonal, linear phase
intMethodStr = {'nearest','bilinear','bicubic'};
level = {1, 2, 3, 4};
n_iterations = length(intMethodStr)*length(waveletStr)*length(level);
alpha = 2;

% Lena Image Load (uint8)
X = imread('lena_std.tif');
img_original = rgb2gray(X);

% Data Holders
rpsnr_down = zeros(1,n_iterations);
rssim_down = zeros(1,n_iterations);
rpsnr_up = zeros(1,n_iterations);
rssim_up = zeros(1,n_iterations);
rniqe = zeros(1,n_iterations);
caseStr = cell(1, n_iterations);
iCount = 0;

% Processing Loop
for ii=1:length(waveletStr)
    for jj = 1:length(level)
    
        % Perform Decomposition per iteration level
        switch level{jj}
            
            % Level 1 DWT ST
            case 1
                for kk=1:length(intMethodStr)
                    
                    iCount = iCount+1;
                    caseStr{iCount} = sprintf('level %d %s - %s',level{jj},waveletStr{ii},intMethodStr{kk});
                    [LL,LH,HL,HH] = dwt2(img_original,waveletStr{ii});
                    
                    jHL = imresize(HL, alpha, intMethodStr{kk});
                    jLH = imresize(LH, alpha, intMethodStr{kk});
                    jHH = imresize(HH, alpha, intMethodStr{kk});
                    
                    img_rec1 = idwt2(img_original,jLH,jHL,jHH,waveletStr{ii},size(img_original).*alpha);
                    img_rec1 = uint8(rescale(img_rec1,0,255));
                    
                    dispBands = [LL LH; HL HH];
                    getDWTSRPlots(img_original, img_rec1, dispBands, waveletStr{ii}, level{jj}, intMethodStr{kk}, alpha);
                    [rpsnr_down(iCount), rpsnr_up(iCount), rssim_down(iCount), rssim_up(iCount), rniqe(iCount)] = ...
                        getDWTSRStats(img_original, img_rec1, intMethodStr{kk}, alpha);
                end % End interpolation method loop
               
                
            % Level 2 DWT ST
            case 2
                for kk=1:length(intMethodStr)
                    
                    iCount = iCount+1;
                    caseStr{iCount} = sprintf('level %d %s - %s',level{jj},waveletStr{ii},intMethodStr{kk});
                    
                    % Level 1
                    [LL,LH,HL,HH] = dwt2(img_original,waveletStr{ii});
                    
                    jHL = imresize(HL, alpha, intMethodStr{kk});
                    jLH = imresize(LH, alpha, intMethodStr{kk});
                    jHH = imresize(HH, alpha, intMethodStr{kk});
                    
                    img_rec1 = idwt2(img_original,jLH,jHL,jHH,waveletStr{ii},size(img_original).*alpha);
                  
                    % Level 2
                    img_rec = imresize(img_rec1,1/alpha,intMethodStr{kk});
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii});

                    jLLHL = imresize(LLHL, alpha, intMethodStr{kk});
                    jLLLH = imresize(LLLH, alpha, intMethodStr{kk});
                    jLLHH = imresize(LLHH, alpha, intMethodStr{kk});

                    img_rec2 = idwt2(img_rec,jLLLH,jLLHL,jLLHH,waveletStr{ii},size(img_original).*alpha);
                    img_rec1 = uint8(rescale(img_rec1,0,255));
                    img_rec2 = uint8(rescale(img_rec2,0,255));
                    
                    dispBands = [LLLL LLLH; LLHL LLHH];
                    getDWTSRPlots(img_original, img_rec2, dispBands, waveletStr{ii}, level{jj}, intMethodStr{kk}, alpha);
                    [rpsnr_down(iCount), rpsnr_up(iCount), rssim_down(iCount), rssim_up(iCount), rniqe(iCount)] = ...
                        getDWTSRStats(img_original, img_rec2, intMethodStr{kk}, alpha);
                end
                
            % Level 3 DWT ST
            case 3
                for kk=1:length(intMethodStr)
                    
                    iCount = iCount+1;
                    caseStr{iCount} = sprintf('level %d %s - %s',level{jj},waveletStr{ii},intMethodStr{kk});
                    
                    % Level 1
                    [LL,LH,HL,HH] = dwt2(img_original,waveletStr{ii});
                    
                    jHL = imresize(HL, alpha, intMethodStr{kk});
                    jLH = imresize(LH, alpha, intMethodStr{kk});
                    jHH = imresize(HH, alpha, intMethodStr{kk});
                    
                    img_rec1 = idwt2(img_original,jLH,jHL,jHH,waveletStr{ii},size(img_original).*alpha);
                  
                    % Level 2
                    img_rec = imresize(img_rec1,1/alpha,intMethodStr{kk});
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii});

                    jLLHL = imresize(LLHL, alpha, intMethodStr{kk});
                    jLLLH = imresize(LLLH, alpha, intMethodStr{kk});
                    jLLHH = imresize(LLHH, alpha, intMethodStr{kk});

                    img_rec2 = idwt2(img_rec,jLLLH,jLLHL,jLLHH,waveletStr{ii},size(img_original).*alpha);
                    
                    % Level 3
                    img_rec = imresize(img_rec2,1/alpha,intMethodStr{kk});
                    [LLLLLL,LLLLLH,LLLLHL,LLLLHH] = dwt2(img_rec,waveletStr{ii});

                    jLLLLHL = imresize(LLLLHL, alpha, intMethodStr{kk});
                    jLLLLLH = imresize(LLLLLH, alpha, intMethodStr{kk});
                    jLLLLHH = imresize(LLLLHH, alpha, intMethodStr{kk});

                    img_rec3 = idwt2(img_rec,jLLLLLH,jLLLLHL,jLLLLHH,waveletStr{ii},size(img_original).*alpha);
                    
                    img_rec1 = uint8(rescale(img_rec1,0,255));
                    img_rec2 = uint8(rescale(img_rec2,0,255));
                    img_rec3 = uint8(rescale(img_rec3,0,255));
                    
                    dispBands = [LLLLLL LLLLLH; LLLLHL LLLLHH];
                    getDWTSRPlots(img_original, img_rec3, dispBands, waveletStr{ii}, level{jj}, intMethodStr{kk}, alpha);
                    [rpsnr_down(iCount), rpsnr_up(iCount), rssim_down(iCount), rssim_up(iCount), rniqe(iCount)] = ...
                        getDWTSRStats(img_original, img_rec3, intMethodStr{kk}, alpha);
                end
                
                % Level 4 DWT ST
            case 4
                for kk=1:length(intMethodStr)
                    
                    iCount = iCount+1;
                    caseStr{iCount} = sprintf('level %d %s - %s',level{jj},waveletStr{ii},intMethodStr{kk});
                    
                    % Level 1
                    [LL,LH,HL,HH] = dwt2(img_original,waveletStr{ii});
                    
                    jHL = imresize(HL, alpha, intMethodStr{kk});
                    jLH = imresize(LH, alpha, intMethodStr{kk});
                    jHH = imresize(HH, alpha, intMethodStr{kk});
                    
                    img_rec1 = idwt2(img_original,jLH,jHL,jHH,waveletStr{ii},size(img_original).*alpha);
                  
                    % Level 2
                    img_rec = imresize(img_rec1,1/alpha,intMethodStr{kk});
                    [LLLL,LLLH,LLHL,LLHH] = dwt2(img_rec,waveletStr{ii});

                    jLLHL = imresize(LLHL, alpha, intMethodStr{kk});
                    jLLLH = imresize(LLLH, alpha, intMethodStr{kk});
                    jLLHH = imresize(LLHH, alpha, intMethodStr{kk});

                    img_rec2 = idwt2(img_rec,jLLLH,jLLHL,jLLHH,waveletStr{ii},size(img_original).*alpha);
                    
                    % Level 3
                    img_rec = imresize(img_rec2,1/alpha,intMethodStr{kk});
                    [LLLLLL,LLLLLH,LLLLHL,LLLLHH] = dwt2(img_rec,waveletStr{ii});

                    jLLLLHL = imresize(LLLLHL, alpha, intMethodStr{kk});
                    jLLLLLH = imresize(LLLLLH, alpha, intMethodStr{kk});
                    jLLLLHH = imresize(LLLLHH, alpha, intMethodStr{kk});

                    img_rec3 = idwt2(img_rec,jLLLLLH,jLLLLHL,jLLLLHH,waveletStr{ii},size(img_original).*alpha);
                    
                    % Level 4
                    img_rec = imresize(img_rec3,1/alpha,intMethodStr{kk});
                    [LLLLLLLL,LLLLLLLH,LLLLLLHL,LLLLLLHH] = dwt2(img_rec,waveletStr{ii});

                    jLLLLLLHL = imresize(LLLLLLHL, alpha, intMethodStr{kk});
                    jLLLLLLLH = imresize(LLLLLLLH, alpha, intMethodStr{kk});
                    jLLLLLLHH = imresize(LLLLLLHH, alpha, intMethodStr{kk});

                    img_rec4 = idwt2(img_rec,jLLLLLLLH,jLLLLLLHL,jLLLLLLHH,waveletStr{ii},size(img_original).*alpha);
                    
                    img_rec1 = uint8(rescale(img_rec1,0,255));
                    img_rec2 = uint8(rescale(img_rec2,0,255));
                    img_rec3 = uint8(rescale(img_rec3,0,255));
                    img_rec4 = uint8(rescale(img_rec4,0,255));
                    
                    dispBands = [LLLLLLLL LLLLLLLH; LLLLLLHL LLLLLLHH];
                    getDWTSRPlots(img_original, img_rec4, dispBands, waveletStr{ii}, level{jj}, intMethodStr{kk}, alpha);
                    [rpsnr_down(iCount), rpsnr_up(iCount), rssim_down(iCount), rssim_up(iCount), rniqe(iCount)] = ...
                        getDWTSRStats(img_original, img_rec4, intMethodStr{kk}, alpha);
                end
                
            otherwise
                error('Unsupported Level Number %d', level{jj});
        end % End level case switch
    end % End level iteration loop
end % End wavelet loop

if ~exist(fullfile('out','stats'),'dir')
    mkdir(fullfile('out','stats'));
end

scrnSize = get(0, 'screensize');

figure; set(gcf, 'Position', scrnSize);
plot((1:1:n_iterations), rpsnr_down); grid on; hold on;
plot((1:1:n_iterations), rpsnr_up);
xlim([1 n_iterations]);
set(gca, 'XTick', 1:n_iterations, 'XTickLabel', caseStr);
xtickangle(90);
legend('PSNR - Downsampled Output','PSNR - Upsampled Input','location','best');
title('Total PSNR'); ylabel('PSNR'); xlabel('Case Number');
fname1 = 'psnrtotal.fig';
savefig(fullfile('out','stats',fname1));

figure; set(gcf, 'Position', scrnSize);
plot((1:1:n_iterations), rssim_down); grid on; hold on;
plot((1:1:n_iterations), rssim_up);
xlim([1 n_iterations]);
set(gca, 'XTick', 1:n_iterations, 'XTickLabel', caseStr);
xtickangle(90);
legend('SSIM - Downsampled Output','SSIM - Upsampled Input','location','best');
title('Total SSIM'); ylabel('SSIM'); xlabel('Case Number');
fname2 = 'ssimtotal.fig';
savefig(fullfile('out','stats',fname2));

figure; set(gcf, 'Position', scrnSize);
plot((1:1:n_iterations), rniqe); grid on; 
xlim([1 n_iterations]);
set(gca, 'XTick', 1:n_iterations, 'XTickLabel', caseStr);
xtickangle(90);
title('NIQE'); ylabel('NIQE'); xlabel('Case Number');
fname3 = 'niqe.fig';
savefig(fullfile('out','stats',fname3));


save(fullfile('out','stats','stats.mat'), 'rpsnr_down', 'rpsnr_up', 'rssim_down', 'rssim_up', 'rniqe');