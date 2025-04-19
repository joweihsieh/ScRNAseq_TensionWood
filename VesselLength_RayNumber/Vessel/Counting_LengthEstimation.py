import numpy as np  # math
import matplotlib.pyplot as plt  # plot graph
from PIL import Image, ImageFilter  # load and filter image
from scipy.ndimage.morphology import binary_erosion  # erode mask
from scipy.ndimage import sobel  # sobel filter for drawing outline
from skimage.measure import label, regionprops  # measuring and label regions
import pandas as pd
import os
import glob

folder_path = '/Users/joweihsieh/Dropbox/YCL/tension_wood/Vessel_size/Vessel_20241230/Raw/'
files = glob.glob(os.path.join(folder_path, '*.tif')) 

valid_regions_list = [] 
for file in files:
    # Load the image
    img = Image.open(file)
    img = img.resize(np.array(img.size) * 3)  # resize image
    img2 = img.filter(ImageFilter.GaussianBlur(1))  # apply Gaussian blur
    img3 = img2.filter(ImageFilter.UnsharpMask(10, 100, 0))  # apply unsharp mask
    
    # Convert image to numpy array
    imgarr = np.array(img3)
    cutoff = 240  # fixed threshold for segmentation
    segment = imgarr[:, :, 2] > cutoff  # create binary segment mask
    segment_erode = binary_erosion(segment, np.ones([1, 1]))  # erode mask
    
    # Label connected components
    labeled = label(segment_erode)
    regions = regionprops(labeled)
    
    # Image dimensions
    img_height, img_width = imgarr.shape[:2]
    
    # Filter valid regions and check if touching the border
    valid_regions = []
    for region in regions:
        if region.eccentricity < 0.98 and region.major_axis_length >= 30:
            min_row, min_col, max_row, max_col = region.bbox
            is_on_border = (
                min_row == 0 or min_col == 0 or
                max_row == img_height or max_col == img_width
            )
            valid_regions.append((region, is_on_border))
    
    # Count valid regions and append to the list
    valid_regions_list.append((file, len(valid_regions)))

    # Extract data for valid regions
    centroids = [region.centroid for region, _ in valid_regions]
    major_axis_lengths = np.array([region.major_axis_length for region, _ in valid_regions]) / 2  # radius
    is_on_border_list = [is_on_border for _, is_on_border in valid_regions]

    # Visualization: Scatter plot of centroids
    plt.figure(figsize=(16, 6))
    plt.imshow(img)

    # Separate centroids based on border status
    centroids_array = np.array(centroids)
    border_indices = np.array(is_on_border_list)
    
    # Plot border cells in gray
    if np.any(border_indices):
        plt.scatter(
            *centroids_array[border_indices].T[::-1], 
            c='gray', 
            s=major_axis_lengths[border_indices] * 6  # scale size
        )

    # Plot non-border cells with color mapping
    if np.any(~border_indices):
        plt.scatter(
            *centroids_array[~border_indices].T[::-1],
            c=major_axis_lengths[~border_indices],
            s=major_axis_lengths[~border_indices] * 6,
            cmap='jet',
            vmin=0, vmax=140
        )
    
    plt.colorbar(label='Major Axis Length (px)')
    plt.title(f"count={len(valid_regions)}")
    plt.tight_layout()
    plt.savefig(f'./Out_series_border_grey/{file.split("/")[-1].split(".")[-2]+"."+str(cutoff)+".optimal.scatter.tif"}')
    
    # Histogram of major axis lengths
    plt.figure(figsize=(10, 8))
    plt.hist(major_axis_lengths, bins=50)
    plt.grid()
    plt.xlabel('Major axis length (px)')
    plt.ylabel('Count')
    plt.savefig(f'./Out_series_border_grey/{file.split("/")[-1].split(".")[-2]+"."+str(cutoff)+".optimal.count.tif"}')

    # Save region data to CSV
    output_df = pd.DataFrame({
        "major_axis_length": major_axis_lengths,
        "is_on_border": is_on_border_list
    })
    output_df.to_csv(
        f'./Out_series_border_grey/{file.split("/")[-1].split(".")[-2]+"."+str(cutoff)+".optimal.count.txt"}',
        sep='\t', index=False
    )
    
    print(file)
    print(valid_regions_list)

# Save overall results to a text file
with open('valid_regions_numbers.txt', 'w') as file:
    for filename, length in valid_regions_list:
        file.write(f'{filename}: {length}\n')
