import argparse
import os
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.plot as plot
import logging
import sys
from datetime import datetime
import matplotlib
import re
import tifffile as tif
matplotlib.use('Agg')

def main(args):
    # Setup logging
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"bigfish_run_all_{args.Date}_{args.Transcript}_{args.Condition}_{current_time}.log"
    log_dir = os.path.join(args.src_dir, "logs", "bigfish_run_all")
    os.makedirs(log_dir, exist_ok=True)
    logging.basicConfig(filename=os.path.join(log_dir, log_filename), 
                            level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"Command: {' '.join(sys.argv)}")
    logging.info(f"Arguments: {args}")

    # Add handler to log messages to console as well
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(console_handler)

    RNA_file_paths = []
    RNA_file_names = []
    
    date_dir = os.path.join(args.src_dir, args.Date)
    if os.path.exists(date_dir) == False:
        exit()
    trans_dir = os.path.join(date_dir, args.Transcript)
    if os.path.exists(trans_dir) == False:
        exit()
    cond_dir = os.path.join(trans_dir, args.Condition)
    if os.path.exists(cond_dir) == False:
        exit()
    stack_dir = os.path.join(cond_dir, "stacks")
    if os.path.exists(stack_dir) == False:
        exit()

    spot_raw_dir = os.path.join(cond_dir,"spots_raw")
    os.makedirs(spot_raw_dir,exist_ok=True)
    spot_raw_plots_dir = os.path.join(cond_dir,"spot_raw_plots")
    os.makedirs(spot_raw_plots_dir,exist_ok=True)

    spot_decomp_dir = os.path.join(cond_dir,"spots_decomp")
    os.makedirs(spot_decomp_dir,exist_ok=True)
    spot_decomp_plots_dir = os.path.join(cond_dir,"spot_decomp_plots")
    os.makedirs(spot_decomp_plots_dir,exist_ok=True)

    for file_name in os.listdir(stack_dir):
        # find files that match the expected format from rename_files.ipynb
        if args.filt_focus:
            pattern = re.compile(args.Date + "~" + args.Transcript + "~" + args.Condition + "~(.*)" + args.rna_chan + "_filtfocus\.tif")
        else:
            pattern = re.compile(args.Date + "~" + args.Transcript + "~" + args.Condition + "~(.*)" + args.rna_chan + "\.tif")
        if pattern.search(file_name):
            RNA_file_names.append(file_name)
            file = os.path.join(stack_dir,file_name)
            RNA_file_paths.append(file)
    
    vox_siz_zyx = (args.vox_z, args.vox_y, args.vox_x)

    for i,file in enumerate(RNA_file_paths):
        im = tif.imread(file)
        spots, threshold = detection.detect_spots(
            images=im,
            threshold=args.threshold, 
            return_threshold=True, 
            voxel_size=vox_siz_zyx,  # in nanometer (one value per dimension zyx)
            spot_radius=tuple([x+50 for x in vox_siz_zyx]))
        logging.info("detected spots for {}".format(RNA_file_names[i]))
        logging.info("\r shape: {0}".format(spots.shape))

        mip = stack.maximum_projection(im)

        raw_plot_file = os.path.join(spot_raw_plots_dir,RNA_file_names[i].replace(".tif","_thresh"+str(args.threshold)+"_raw_spots"))
        plot.plot_detection(mip, spots, contrast=True,path_output=raw_plot_file,show=False)

        raw_spot_file = os.path.join(spot_raw_dir,RNA_file_names[i].replace(".tif","_thresh"+str(args.threshold)+"_raw_spots.csv"))
        stack.save_data_to_csv(spots,raw_spot_file)

        # now run decomposition
        spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
            image=im, 
            spots=spots, 
            voxel_size=vox_siz_zyx,  # in nanometer (one value per dimension zyx)
            spot_radius=tuple([x+50 for x in vox_siz_zyx]),
            alpha=0.7,  # alpha impacts the number of spots per candidate region
            beta=1,  # beta impacts the number of candidate regions to decompose
            gamma=5)  # gamma the filtering step to denoise the image
        logging.info("detected spots after decomposition for {}".format(RNA_file_names[i]))
        logging.info("\r shape: {0}".format(spots_post_decomposition.shape))

        decomp_plot_file = os.path.join(spot_decomp_plots_dir,RNA_file_names[i].replace(".tif","_thresh"+str(args.threshold)+"_decomp_spots"))
        plot.plot_detection(mip, spots_post_decomposition, contrast=True,path_output=decomp_plot_file,show=False)

        decomp_spot_file = os.path.join(spot_decomp_dir,RNA_file_names[i].replace(".tif","_thresh"+str(args.threshold)+"_decomp_spots.csv"))
        stack.save_data_to_csv(spots_post_decomposition,decomp_spot_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run detection on all images in the src_dir.')
    parser.add_argument('Date', type=str, help='Date of the experiment.')
    parser.add_argument('Transcript', type=str, help='Transcript name.')
    parser.add_argument('Condition', type=str, help='Experimental condition.')
    parser.add_argument('src_dir', type=str, help='Source directory containing the data.')
    parser.add_argument('rna_chan', type=str, help='RNA channel name.')
    parser.add_argument('vox_x', type=int, help='Voxel size in X dimension.')
    parser.add_argument('vox_y', type=int, help='Voxel size in Y dimension.')
    parser.add_argument('vox_z', type=int, help='Voxel size in Z dimension.')
    parser.add_argument('threshold', type=int, help='Threshold value for detection.')
    parser.add_argument('--filt_focus', action='store_true', default=False, help='Whether to use filtered images.')
    args = parser.parse_args()

    main(args)
