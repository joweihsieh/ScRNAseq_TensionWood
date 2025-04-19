import numpy as np # math
import matplotlib.pyplot as plt # plot graph
from PIL import Image, ImageFilter # load and filter image
from scipy.ndimage.morphology import binary_erosion # erode mask
from scipy.ndimage import sobel # sobel filter for drawing outline
from skimage.measure import label, regionprops # measuring and label regions
import pandas as pd
import os
import glob


#files = [
#	'./Raw/peptide_8_4th_N.tiff',
#	'./Raw/peptide8_4th_1uM.tiff'
#]

folder_path = '/Users/joweihsieh/Dropbox/YCL/tension_wood/ray_number/20250109_ray quantification/Raw/'
files = glob.glob(os.path.join(folder_path, '*.tif')) 

#file = "/Users/joweihsieh/Dropbox/YCL/tension_wood/ray_number/Raw/Ptr_vertical_bio3_20X-3_colored.jpg"

valid_regions_list = [] 
for file in files:
    number = files.index(file)
    cutoff = 240 # fixed !!! 

    img = Image.open(file)
    img = img.resize(np.array(img.size)*3)
    img2 = img.filter(ImageFilter.GaussianBlur(1))
    img3 = img2.filter(ImageFilter.UnsharpMask(10, 100, 0))


    imgarr = np.array(img3)
    segment = imgarr[:, :, 2]> cutoff
    segment_erode = binary_erosion(segment, np.ones([1, 1]))


    labeled = label(segment_erode)
    regions = regionprops(labeled)
    eccentricity = [region.eccentricity for region in regions] 
    valid_regions = [region for region in regions if region.eccentricity < 1 and region.major_axis_length >= 20] 
    #valid_regions_list.append(len(valid_regions))
    valid_regions_list.append((file, len(valid_regions)))

    centroid = [region.centroid for region in valid_regions]
    major_axis_length = np.array([region.major_axis_length for region in valid_regions])
    minor_axis_length = np.array([region.minor_axis_length for region in valid_regions])

    radius_major = major_axis_length / 2  
    radius_minor = minor_axis_length / 2  


    plt.figure(figsize=(16, 6))
    plt.imshow(img)
    plt.scatter(*np.array(centroid).T[::-1], c=radius_major, s=radius_major*6, cmap='jet', vmin=0, vmax=140)
    plt.colorbar()
    plt.title(f"count={len(valid_regions)}")
    plt.tight_layout()

    plt.savefig(f'./Out_series/{file.split("/")[-1].split(".")[-2]+"."+str(cutoff)+".optimal.scatter.tif"}')
    print("max cell number")
    print(file)



    plt.figure(figsize=(10, 8))
    plt.hist(major_axis_length, bins=50)

    plt.grid()
    plt.xlabel('major axis length (px)')
    plt.ylabel('count')
    plt.savefig(f'./Out_series/{file.split("/")[-1].split(".")[-2]+"."+str(cutoff)+".optimal.count.tif"}')

    my_df = pd.DataFrame(major_axis_length)
    my_df.to_csv(f'./Out_series/{file.split("/")[-1].split(".")[-2]+"."+str(cutoff)+".optimal.count.txt"}',sep='\t',index=False,header=False)

    print(file)
    print(valid_regions_list)

with open('valid_regions_numbers.txt', 'w') as file:
    for length in valid_regions_list:
        file.write(f'{length}\n')


with open('valid_regions_lengths.txt', 'w') as file:
    for filename, length in valid_regions_list:
        file.write(f'{filename}: {length}\n')


        



