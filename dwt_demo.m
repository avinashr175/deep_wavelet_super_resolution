clear all; close all;

load woman;
img_original = X;

figure; imagesc(img_original); colormap(gray);
title(sprintf('%d x %d Input Image',size(img_original,1),size(img_original,2)));

[LL,LH,HL,HH] = dwt2(img_original,'db1');

img_disp = [LL HL; LH HH];
figure; imagesc(img_disp); colormap(gray);
title('Daubechies Level 1 Sub-band Decomposition');

[LLLL,LLLH,LLHL,LLHH] = dwt2(LL,'db1');

ll_disp = [LLLL LLHL; LLLH LLHH];
img_disp = [ll_disp HL; LH HH];
figure; imagesc(img_disp); colormap(gray);
title('Daubechies Level 2 Sub-band Decomposition');

jHL = imresize(HL, 2, 'bicubic');
jLH = imresize(LH, 2, 'bicubic');
jHH = imresize(HH, 2, 'bicubic');

img_rec = idwt2(img_original,jLH,jHL,jHH,'db1',size(img_original).*2);
figure; imagesc(img_rec); colormap(gray);
title(sprintf('%d x %d 1 Level Reconstructed Image',size(img_rec,1),size(img_rec,2)));

jLLHL = imresize(LLHL, 2*2, 'bicubic');
jLLLH = imresize(LLLH, 2*2, 'bicubic');
jLLHH = imresize(LLHH, 2*2, 'bicubic');

img_rec = idwt2(img_original,jLLHL,jLLHL,jLLHH,'db1',size(img_original).*2);
figure; imagesc(img_rec); colormap(gray);
title(sprintf('%d x %d 2 Level Reconstructed Image',size(img_rec,1),size(img_rec,2)));

