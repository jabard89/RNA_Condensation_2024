import os
import tifffile
import numpy as np
import argparse

def combine_masks(directory):
    # Use os.path.basename() to get the terminal directory's name
    folder_name = os.path.basename(directory)
    combined_mask_name = os.path.join(directory, folder_name + "_masks.tiff")
    
    # Your existing logic to combine the masks
    file_list = sorted([f for f in os.listdir(directory) if '_Cell_' in f])
    combined_mask = None

    for file_name in file_list:
        cell_id = int(file_name.split('_Cell_')[-1].split('.tiff')[0])
        mask_path = os.path.join(directory, file_name)
        current_mask = tifffile.imread(mask_path)

        if combined_mask is None:
            combined_mask = np.zeros_like(current_mask)

        combined_mask[current_mask > 0] = cell_id

    tifffile.imwrite(combined_mask_name, combined_mask)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine individual cell masks into a single mask.')
    parser.add_argument('directory', type=str, help='Directory containing the individual cell masks.')
    args = parser.parse_args()
    combine_masks(args.directory)
