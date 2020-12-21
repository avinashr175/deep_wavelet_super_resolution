import argparse
import glob
import h5py
import numpy as np
import PIL.Image as pil_image
from utils import convert_rgb_to_y

h5_file = h5py.File('data.h5', 'w')

images_dir = 'DIV2K_train_LR_bicubic/hr_images'
scale = 2
patch_size = 50
stride = 30

lr_patches = []
hr_patches = []

count=0

for image_path in sorted(glob.glob('{}/*'.format(images_dir))):
    count+=1
    if(count==100):
        break
    print(image_path)
    hr = pil_image.open(image_path).convert('RGB')
    hr_width = (hr.width // scale) * scale
    hr_height = (hr.height // scale) * scale
    hr = hr.resize((hr_width, hr_height), resample=pil_image.BICUBIC)
    lr = hr.resize((hr_width // scale, hr_height // scale), resample=pil_image.BICUBIC)
    lr = lr.resize((lr.width * scale, lr.height * scale), resample=pil_image.BICUBIC)
    hr = np.array(hr).astype(np.float32)
    lr = np.array(lr).astype(np.float32)
    hr = convert_rgb_to_y(hr)
    lr = convert_rgb_to_y(lr)

    for i in range(0, lr.shape[0] - patch_size + 1, stride):
        for j in range(0, lr.shape[1] - patch_size + 1, stride):
            lr_patches.append(lr[i:i + patch_size, j:j + patch_size])
            hr_patches.append(hr[i:i + patch_size, j:j + patch_size])

lr_patches = np.array(lr_patches)
hr_patches = np.array(hr_patches)
print("hi")
h5_file.create_dataset('lr', data=lr_patches)
h5_file.create_dataset('hr', data=hr_patches)

h5_file.close()