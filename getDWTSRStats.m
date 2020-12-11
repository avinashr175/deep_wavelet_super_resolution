function [psnr_down, psnr_up, ssim_down, ssim_up, inique] = getDWTSRStats(og, rec, intStr, alpha)

rec_dwn = imresize(rec, 1/alpha, intStr);
psnr_down = psnr(og, rec_dwn);
ssim_down = ssim(og, rec_dwn);

og_up = imresize(og,alpha,intStr);
psnr_up = psnr(og_up, rec);
ssim_up = ssim(og_up, rec);

inique = niqe(rec);

end