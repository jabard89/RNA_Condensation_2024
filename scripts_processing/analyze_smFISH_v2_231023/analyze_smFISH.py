
import os
import tifffile as tif
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import numpy as np
import pandas as pd
import re
import warnings
from scipy.ndimage import binary_erosion
import argparse
import logging
from datetime import datetime
import sys
matplotlib.use('Agg')

def setup_logging(log_dir, log_name, args):
    # Configure the root logger
    logging.basicConfig(filename=os.path.join(log_dir, log_name), 
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    # Get the root logger
    logger = logging.getLogger()
    
    # Log command and arguments
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info(f"Arguments: {args}")

    # Add handler to log messages to console as well
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(console_handler)
    
    return logger

def match_point_to_cell_2D(mask,location,radius=2):
    y, x = location
    rows, cols = mask.shape
    
    # Check if location is within the bounds of the array
    if y < 0 or y >= rows or x < 0 or x >= cols:
        warnings.warn(f"Location {location} is outside the bounds of the array")
        return None
    
    value = mask[y, x]
    # Check if location is not part of an object (has a zero value)
    if value == 0:
        return None
    # the binary_erosion function shrinks the true objects by basically rolling a ball around the edges of the object
    eroded_mask = binary_erosion(mask == value, structure=np.ones((radius+1, radius+1)))
    if eroded_mask[y,x] == False:
        return None
    # If all conditions are met, return the value of the object
    return value

def load_cell_norm_2D(image, mask, cell_id):
    # load an image and extract the region of interest, intensities normalized to the mean intensity of the cell
    # first get the location of the cell
    # set values outside of the mask to nan, normalize values inside of the cell to the mean intensity of the cell
    # uses 2D mask and 2d image
    if mask.shape != image.shape:
        raise ValueError("Image and mask must be the same size")
    image = image.astype(np.float32)
    image[mask != cell_id] = np.nan
    if np.isnan(image).all():
        return np.zeros(image.shape)
    image = image / np.nanmean(image)
    image[np.isnan(image)] = 0
    return(image)

def extract_square(array,center,radius=1):
    y, x = center
    rows, cols = array.shape
    
    # Check if location is within the bounds of the array
    if y < radius or y >= (rows-radius) or x < radius or x >= (cols-radius):
        warnings.warn(f"Location {center} is outside the bounds of the array")
        return None
    return(array[y-radius:y+radius,x-radius:x+radius])

def plot_rna_slice(img,points,output_file,buffer=1):
    # Plot image
    fig, ax = plt.subplots()
    ax.imshow(img, cmap='viridis')  # Change colormap to 'viridis'
    side_length = 2*buffer + 1
    # Draw squares centered on points
    for y, x in points:
        # Determine lower left corner of square
        bottom_left_x = x - buffer
        bottom_left_y = y - buffer

        # Draw square (already set to red in the previous code)
        rect = matplotlib.patches.Rectangle((bottom_left_x, bottom_left_y), side_length, side_length,
                                linewidth=0.25, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    # Save as PNG
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

    return(None)

def get_valid_indices(mask, object_id, radius):
    def inShape(loc,shape,radius):
        y = loc[0]
        x = loc[1]
        if (x-radius) < 0 or (x+radius) >= shape[1]:
            return(False)
        if (y-radius) < 0 or (y+radius) >= shape[0]:
            return(False)
        return(True)
    # Erode the mask to include only the area at least `buffer` distance away from the border.
    eroded_mask = binary_erosion(mask == object_id,structure=np.ones((radius+1, radius+1)))
    
    # Extract the y, x coordinates of the valid points (where eroded_mask is True)
    valid_y, valid_x = np.where(eroded_mask)
    locations = [loc for loc in zip(valid_y, valid_x) if inShape(loc,mask.shape,radius)]
    
    return locations

def get_random_points_from_valid_indices(valid_indices, n):
    # Randomly select n points from the list of valid indices.
    # Note: If n > len(valid_indices), you should either reduce n or enable replacement in np.random.choice
    selected_indices = np.random.choice(len(valid_indices), min(n, len(valid_indices)), replace=True)
    
    return [valid_indices[i] for i in selected_indices]

def create_rna_dataframe(fov_dict):
    rna_df = pd.DataFrame(fov_dict['RNA_spots']['RNA_id'],columns=['RNA_id'])
    rna_df['Cell_id'] = fov_dict['RNA_spots']['Cell_id']
    (loc_y, loc_x) = [x[0] for x in fov_dict['RNA_spots']['Location']], [x[1] for x in fov_dict['RNA_spots']['Location']]
    rna_df['RNA_loc_y'] = loc_y
    rna_df['RNA_loc_x'] = loc_x
    rna_df['Zidx'] = fov_dict['RNA_spots']['Zidx']
    rna_df['Transcript'] = fov_dict['trans']
    rna_df['Condition'] = fov_dict['cond']
    rna_df['Date'] = fov_dict['date']
    rna_df['Image_id'] = fov_dict['image_id']
    
    mean_PAB1_list = []
    mean_RNA_list = []
    for index, row in rna_df.iterrows():
        if row['Cell_id'] is not None and not np.isnan(row['Cell_id']):
            cell_id = row['Cell_id']
            Zidx = row['Zidx']
            RNA_id = row['RNA_id']
            mean_PAB1_list.append(fov_dict['cells'][cell_id][Zidx]['RNA'][RNA_id]['mean_PAB1'])
            mean_RNA_list.append(fov_dict['cells'][cell_id][Zidx]['RNA'][RNA_id]['mean_RNA'])
        else:
            mean_PAB1_list.append(np.nan)
            mean_RNA_list.append(np.nan)
    rna_df['RNA_mean_PAB1'] = mean_PAB1_list
    rna_df['RNA_mean_RNA'] = mean_RNA_list
    return(rna_df)

# extract long tidy dataframe of RNA_mean_PAB1 data
def create_random_box_df(fov_dict):
    df_rows = []
    Date = fov_dict['date']
    Transcript = fov_dict['trans']
    Condition = fov_dict['cond']
    Image_id = fov_dict['image_id']
    for cell_id,cell in fov_dict['cells'].items():
        for Zidx,cell_slice in cell.items():
            for i,location in enumerate(cell_slice['random_points']):
                loc_y = location[0]
                loc_x = location[1]
                random_mean_PAB1 = cell_slice['random_points_mean_PAB1'][i]
                random_mean_RNA = cell_slice['random_points_mean_RNA'][i]
                row = {'Date':Date,'Transcript':Transcript,'Condition':Condition,'Image_id':Image_id,'Cell_id':cell_id,
                    'Zidx':Zidx,'random_loc_y':loc_y,'random_loc_x':loc_x,'random_mean_PAB1':random_mean_PAB1,'random_mean_RNA':random_mean_RNA}
                df_rows.append(row)
    return(pd.DataFrame(df_rows))

def main(Date, Transcript, Condition, src_dir, RNA_chan, PAB1_chan, threshold, nRandomPoints, RNA_buffer, filt_focus, args):
    # Convert script content here
    
    date_dir = os.path.join(src_dir, Date)
    if os.path.exists(date_dir) == False:
        exit()
    trans_dir = os.path.join(date_dir, Transcript)
    if os.path.exists(trans_dir) == False:
        exit()
    cond_dir = os.path.join(trans_dir, Condition)
    if os.path.exists(cond_dir) == False:
        exit()
    stack_dir = os.path.join(cond_dir, "stacks")
    if os.path.exists(stack_dir) == False:
        exit()
    mask_dir = os.path.join(src_dir,"combined_masks")
    if os.path.exists(mask_dir) == False:
        exit()
    
    if filt_focus:
        pattern = re.compile(f"{Date}~{Transcript}~{Condition}~(?P<FOV>.+?)~{PAB1_chan}_filtfocus\.tif$")
    else:
        pattern = re.compile(f"{Date}~{Transcript}~{Condition}~(?P<FOV>.+?)~{PAB1_chan}\.tif$")

    # load FOVS
    fovs = {}
    all_rna_df = []
    all_random_df = []

    for file in os.listdir(stack_dir):
        match = pattern.search(file)
        if match:
            image_id = match.group('FOV')
            image_path = os.path.join(stack_dir,file)
            fov_name = '~'.join([Date,Transcript,Condition,image_id])
            if filt_focus:
                mask_path = os.path.join(mask_dir,fov_name+f'~{PAB1_chan}_filtfocus_mip_masks.tiff')
            else:
                mask_path = os.path.join(mask_dir,fov_name+f'~{PAB1_chan}_mip_masks.tiff')
            if os.path.isfile(mask_path):
                decom_spot_file = os.path.join(cond_dir,'spots_decomp',file.replace(PAB1_chan,RNA_chan).replace('.tif',f'_thresh{threshold}_decomp_spots.csv'))
                raw_spot_file = os.path.join(cond_dir,'spots_raw',file.replace(PAB1_chan,RNA_chan).replace('.tif',f'_thresh{threshold}_raw_spots.csv'))
                if os.path.isfile(decom_spot_file):
                    fovs[fov_name] = {"Name":fov_name,"date":Date,"trans":Transcript,"cond":Condition,"image_id":image_id,
                                      "stack_paths":{PAB1_chan:image_path,RNA_chan:image_path.replace(PAB1_chan,RNA_chan)},
                                      "mask_path":mask_path,
                                      "decomp_spot_file":decom_spot_file,"raw_spot_file":raw_spot_file}
    if len(fovs) == 0:
        exit()
    # Set up logging
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = os.path.join(src_dir, "logs", "analyze_smFISH")
    os.makedirs(log_dir, exist_ok=True)
    log_name = f"analyze_smFISH_{Date}_{Transcript}_{Condition}_{current_time}.log"
    logger = setup_logging(log_dir, log_name, args)
    logger.info("Script started")

    # Now lets load the cells and RNAs
    for fov,fov_dict in fovs.items():
        stack_PAB1 = tif.imread(fov_dict['stack_paths'][PAB1_chan])
        stack_RNA = tif.imread(fov_dict['stack_paths'][RNA_chan])
        decom_spot_file = fov_dict['decomp_spot_file']
        raw_spot_file = fov_dict['raw_spot_file']

        # first load the cells
        mask = tif.imread(fov_dict['mask_path'])
        nCells = mask.max()
        nSlice = stack_PAB1.shape[0]
        norm_cells_PAB1 = [ [] for _ in range(nSlice) ]
        norm_cells_RNA = [ [] for _ in range(nSlice) ]
        fov_dict['cells'] = {}
        for i in range(nCells):
            cell_id = i+1
            fov_dict['cells'][cell_id] = {}
            for Zidx in range(nSlice):
                fov_dict['cells'][cell_id][Zidx] = {'nRNA':0,'RNA':{}}
                # preload the normalized cells
                norm_cell_PAB1 = load_cell_norm_2D(stack_PAB1[Zidx,:,:], mask, cell_id)
                norm_cells_PAB1[Zidx].append(norm_cell_PAB1)
                norm_cell_RNA = load_cell_norm_2D(stack_RNA[Zidx,:,:], mask, cell_id)
                norm_cells_RNA[Zidx].append(norm_cell_RNA)
        # combine all images in a stack into one image
        norm_cells_PAB1 = [np.sum(slice,axis=0) for slice in norm_cells_PAB1]
        norm_cells_RNA = [np.sum(slice,axis=0) for slice in norm_cells_RNA]

        # now load the RNA
        logger.info(f"Loading RNA for {raw_spot_file}")
        fov_dict["spots_raw"] = pd.read_csv(raw_spot_file,names=["Z","Y","X"],sep=";").to_numpy()
        fov_dict["spots_decomp"] = pd.read_csv(decom_spot_file,names=["Z","Y","X"],sep=";").to_numpy()
        
        # now associate RNA with cells
        fov_dict["RNA_spots"] = {"RNA_id":[],"Cell_id":[],"Location":[],"Zidx":[]}
        orphan_rna_i = 1
        for row in fov_dict["spots_raw"]:
            Zidx = int(row[0])
            location = (int(row[1]),int(row[2])) # y,x
            cell_id = match_point_to_cell_2D(mask,location,radius=3)
            if cell_id:
                # extract the Pab1 signal in a square around the spot
                # first load the image
                cell_slice = fov_dict['cells'][cell_id][Zidx]
                cell_slice['nRNA'] += 1
                rna_id = f"Cell{cell_id}_Z{Zidx}_RNA{cell_slice['nRNA']}"
                cell_slice['RNA'][rna_id] = {}
                
                rna_box_PAB1 = extract_square(norm_cells_PAB1[Zidx],location,radius=RNA_buffer)
                if rna_box_PAB1 is None:
                    cell_slice['RNA'][rna_id]['mean_PAB1'] = np.nan
                else:
                    rna_box_PAB1[rna_box_PAB1 == 0] = np.nan
                    cell_slice['RNA'][rna_id]['mean_PAB1'] = np.nanmean(rna_box_PAB1)
                
                rna_box_RNA = extract_square(norm_cells_RNA[Zidx],location,radius=RNA_buffer)
                if rna_box_RNA is None:
                    cell_slice['RNA'][rna_id]['mean_RNA'] = np.nan
                else:
                    rna_box_RNA[rna_box_RNA == 0] = np.nan
                    cell_slice['RNA'][rna_id]['mean_RNA'] = np.nanmean(rna_box_RNA)
            else:
                cell_id = None
                rna_id = f"Orphan{orphan_rna_i}"
                orphan_rna_i += 1
            fov_dict["RNA_spots"]["Cell_id"].append(cell_id)
            fov_dict["RNA_spots"]["Location"].append(location)
            fov_dict["RNA_spots"]["Zidx"].append(Zidx)
            fov_dict['RNA_spots']['RNA_id'].append(rna_id)
        
        # create a dictionary of the fov to output
        rna_df = create_rna_dataframe(fov_dict)
        all_rna_df.append(rna_df)

        # now write the overlayed RNAs to an image
        # make a dataframe of the RNA
        RNA_df = pd.DataFrame(fov_dict["RNA_spots"])
        RNA_outline_dir = os.path.join(cond_dir,'RNA_outlines',fov)
        os.makedirs(RNA_outline_dir,exist_ok=True)

        for Zidx in range(len(norm_cells_PAB1)):
            # first filter filter out any orphan RNAs or RNAs that were filtered for being too close to the edge
            RNA_df_Zidx = RNA_df[(RNA_df['Zidx'] == Zidx)]
            RNA_df_Zidx = RNA_df_Zidx[RNA_df_Zidx['Cell_id'].notnull()]
            RNA_locations = RNA_df_Zidx['Location'].tolist()
            out_file_PAB1 = os.path.join(RNA_outline_dir,f"{fov}_Z{Zidx}_RNA_norm{PAB1_chan}.png")
            out_file_RNA = os.path.join(RNA_outline_dir,f"{fov}_Z{Zidx}_RNA_norm{RNA_chan}.png")
            if np.max(norm_cells_PAB1[Zidx]) == 0: # this would happen if there are no cells in the mask slice
                plot_rna_slice(np.zeros(mask.shape),RNA_locations,out_file_PAB1,buffer=RNA_buffer)
                plot_rna_slice(np.zeros(mask.shape),RNA_locations,out_file_RNA,buffer=RNA_buffer)
                continue
            plot_rna_slice(norm_cells_PAB1[Zidx],RNA_locations,out_file_PAB1,buffer=RNA_buffer)
            plot_rna_slice(norm_cells_RNA[Zidx],RNA_locations,out_file_RNA,buffer=RNA_buffer)


        # finally lets calculate random points in the cells
        for cell_id,cell in fov_dict['cells'].items():
            for Zidx,cell_slice in cell.items():
                valid_indices = get_valid_indices(mask,object_id=cell_id,radius=3)
                points_needed = nRandomPoints
                random_points = get_random_points_from_valid_indices(valid_indices, points_needed)
                cell_slice['random_points'] = random_points
                random_boxes_PAB1 = [extract_square(norm_cells_PAB1[Zidx],location,radius=RNA_buffer) for location in random_points]
                random_boxes_RNA = [extract_square(norm_cells_RNA[Zidx],location,radius=RNA_buffer) for location in random_points]
                cell_slice['random_points_mean_PAB1'] = [np.nanmean(box) for box in random_boxes_PAB1]
                cell_slice['random_points_mean_RNA'] = [np.nanmean(box) for box in random_boxes_RNA]

        random_df = create_random_box_df(fov_dict)
        all_random_df.append(random_df)

    # now export the dataframes    
    all_rna_df = pd.concat(all_rna_df)
    all_random_df = pd.concat(all_random_df)
    all_rna_df.to_csv(os.path.join(cond_dir,f"{Date}~{Transcript}~{Condition}_rna_df.tsv.gz"),sep="\t",index=False,compression="gzip")
    all_random_df.to_csv(os.path.join(cond_dir,f"{Date}~{Transcript}~{Condition}_random_df.tsv.gz"),sep="\t",index=False,compression="gzip")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze smFISH data.')
    parser.add_argument('Date', type=str, help='Date of the experiment.')
    parser.add_argument('Transcript', type=str, help='Transcript name.')
    parser.add_argument('Condition', type=str, help='Experimental condition.')
    parser.add_argument('src_dir', type=str, help='Source directory containing the data.')
    parser.add_argument('rna_chan', type=str, help='RNA channel name.')
    parser.add_argument('pab1_chan', type=str, help='PAB1 channel name.')
    parser.add_argument('threshold', type=int, help='Threshold value for detection.')
    parser.add_argument('nRandomPoints', type=int, help='Number of random points per cell to generate.')
    parser.add_argument('RNA_buffer', type=int, help='Buffer around spot for ROI.')
    parser.add_argument('--filt_focus', action='store_true', help='Whether to use filtered images.')
    args = parser.parse_args()

    main(args.Date, args.Transcript, args.Condition, args.src_dir, args.rna_chan, args.pab1_chan, args.threshold, args.nRandomPoints, args.RNA_buffer, args.filt_focus, args)
