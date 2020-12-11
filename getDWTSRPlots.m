function getDWTSRPlots(og, rec, dispBands, wStr, level, intStr, alpha)

if ~exist('out','dir')
    mkdir('out');
end

switch level
    case 1
        % Display Parameters
        scrnSize = get(0, 'screensize');
        dispSize = [1 1 scrnSize(3) scrnSize(4)/2];

        % Save Raw Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,3,1); imagesc(og); title('Original Image');
        subplot(1,3,2); imagesc(dispBands);
        title(sprintf('%s, %s Interpolation\nDWT Sub-Band Decomposition Bands',...
                        wStr, intStr));
        subplot(1,3,3); imagesc(rec); 
        title('Level 1 Reconstruction'); %% fix!
        fname1 = ['recout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname1));
        
        % Save Downsampled Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og); title('Original Image');
        rec_dwn = imresize(rec, 1/alpha, intStr);
        subplot(1,2,2); imagesc(rec_dwn); 
        title(sprintf('%s, %s Interpolation\nDownsampled Level 1 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og, rec_dwn)));
        fname2 = ['dwnrecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname2));
        
        % Save Upsampled Input
        og_up = imresize(og,alpha,intStr);
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og_up); title('Upsampled Original Image');
        subplot(1,2,2); imagesc(rec); 
        title(sprintf('%s, %s Interpolation\nLevel 1 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og_up, rec)));
        fname3 = ['uprecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname3));
    
    case 2
        % Display Parameters
        scrnSize = get(0, 'screensize');
        dispSize = [1 1 scrnSize(3) scrnSize(4)/2];

        % Save Raw Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,3,1); imagesc(og); title('Original Image');
        subplot(1,3,2); imagesc(dispBands);
        title(sprintf('%s, %s Interpolation\nDWT Sub-Band Decomposition Bands',...
                        wStr, intStr));
        subplot(1,3,3); imagesc(rec); 
        title('Level 2 Reconstruction'); %% fix!
        fname1 = ['recout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname1));
        
        % Save Downsampled Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og); title('Original Image');
        rec_dwn = imresize(rec, 1/alpha, intStr);
        subplot(1,2,2); imagesc(rec_dwn); 
        title(sprintf('%s, %s Interpolation\nDownsampled Level 2 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og, rec_dwn)));
        fname2 = ['dwnrecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname2));
        
        % Save Upsampled Input
        og_up = imresize(og,alpha,intStr);
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og_up); title('Upsampled Original Image');
        subplot(1,2,2); imagesc(rec); 
        title(sprintf('%s, %s Interpolation\nLevel 2 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og_up, rec)));
        fname3 = ['uprecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname3));
        
    case 3
        
        % Display Parameters
        scrnSize = get(0, 'screensize');
        dispSize = [1 1 scrnSize(3) scrnSize(4)/2];

        % Save Raw Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,3,1); imagesc(og); title('Original Image');
        subplot(1,3,2); imagesc(dispBands);
        title(sprintf('%s, %s Interpolation\nDWT Sub-Band Decomposition Bands',...
                        wStr, intStr));
        subplot(1,3,3); imagesc(rec); 
        title('Level 3 Reconstruction'); %% fix!
        fname1 = ['recout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname1));
        
        % Save Downsampled Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og); title('Original Image');
        rec_dwn = imresize(rec, 1/alpha, intStr);
        subplot(1,2,2); imagesc(rec_dwn); 
        title(sprintf('%s, %s Interpolation\nDownsampled Level 3 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og, rec_dwn)));
        fname2 = ['dwnrecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname2));
        
        % Save Upsampled Input
        og_up = imresize(og,alpha,intStr);
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og_up); title('Upsampled Original Image');
        subplot(1,2,2); imagesc(rec); 
        title(sprintf('%s, %s Interpolation\nLevel 3 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og_up, rec)));
        fname3 = ['uprecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname3));
        
    case 4
        
        % Display Parameters
        scrnSize = get(0, 'screensize');
        dispSize = [1 1 scrnSize(3) scrnSize(4)/2];

        % Save Raw Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,3,1); imagesc(og); title('Original Image');
        subplot(1,3,2); imagesc(dispBands);
        title(sprintf('%s, %s Interpolation\nDWT Sub-Band Decomposition Bands',...
                        wStr, intStr));
        subplot(1,3,3); imagesc(rec); 
        title('Level 4 Reconstruction'); %% fix!
        fname1 = ['recout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname1));
        
        % Save Downsampled Output
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og); title('Original Image');
        rec_dwn = imresize(rec, 1/alpha, intStr);
        subplot(1,2,2); imagesc(rec_dwn); 
        title(sprintf('%s, %s Interpolation\nDownsampled Level 4 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og, rec_dwn)));
        fname2 = ['dwnrecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname2));
        
        % Save Upsampled Input
        og_up = imresize(og,alpha,intStr);
        figure; set(gcf, 'Position', dispSize);
        subplot(1,2,1); imagesc(og_up); title('Upsampled Original Image');
        subplot(1,2,2); imagesc(rec); 
        title(sprintf('%s, %s Interpolation\nLevel 4 Reconstruction\nPSNR = %d',...
                        wStr, intStr, psnr(og_up, rec)));
        fname3 = ['uprecout_' wStr '_' intStr '_' num2str(level) '.fig'];
        savefig(fullfile('out',fname3));
    
    otherwise
        error(sprintf('Unsupported Level Number %d', level));
end

close all;
end