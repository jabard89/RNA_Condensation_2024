import argparse
import os
import random
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.plot as plot
import logging
import sys
from datetime import datetime
import matplotlib
import re
matplotlib.use('Agg')

def main(args):
    # Setup logging
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"bigfish_test_thresh_{args.Date}_{args.Transcript}_{args.Condition}_{current_time}.log"
    log_dir = os.path.join(args.src_dir, "logs", "bigfish_test_thresh")
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
            # Create "test_threshold" directory
            test_threshold_dir = os.path.join(args.src_dir,args.Date,args.Transcript,args.Condition,"test_threshold")
            os.makedirs(test_threshold_dir, exist_ok=True)
    
    if len(RNA_file_paths) == 0:
        print("No files found.")
        exit()

    
    # Select random subset
    num_samples = min(6, len(RNA_file_paths))
    random_i = random.sample(range(len(RNA_file_paths)), num_samples)
    test_files = [RNA_file_paths[i] for i in random_i]
    test_names = [RNA_file_names[i] for i in random_i]
    test_images = [stack.read_image(file) for file in test_files]
    
    vox_siz_zyx = (args.vox_z, args.vox_y, args.vox_x)

    for i,im in enumerate(test_images):
        elbow_file = os.path.join(test_threshold_dir, test_names[i].replace(".tif", "_elbow"))
        plot.plot_elbow(
            title=test_names[i],
            images=im,
            voxel_size=vox_siz_zyx, 
            spot_radius=tuple([x+50 for x in vox_siz_zyx]),
            path_output=elbow_file,
            show=False
        )
    
    for i, im in enumerate(test_images):
        spots, threshold = detection.detect_spots(
            images=im,
            threshold=args.threshold, 
            return_threshold=True, 
            voxel_size=vox_siz_zyx,
            spot_radius=tuple([x + 50 for x in vox_siz_zyx]))
        spot_file = os.path.join(test_threshold_dir,RNA_file_names[i].replace(".tif",f"_thresh{str(args.threshold)}_spots.csv"))
        stack.save_data_to_csv(spots,spot_file)
        logging.info(f"Detected spots for test image {i}")
        logging.info(f"Shape: {spots.shape}")
        logging.info(f"Saved spots to: {spot_file}")

        detection_file = os.path.join(test_threshold_dir, test_names[i].replace(".tif", f"_thresh{str(args.threshold)}_detected"))
        # Save outputs in the "test_threshold" directory
        mip = stack.maximum_projection(im)
        plot.plot_detection(mip, spots, contrast=True,
            path_output=detection_file,
            show=False)
        logging.info(f"Saved output to {detection_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run detection on a random subset of images.')
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
