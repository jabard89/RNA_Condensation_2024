import argparse
import os
import re
import tifffile
import bigfish.stack as stack
import logging
import sys
from datetime import datetime
import matplotlib
matplotlib.use('Agg')

def pick_most_infocus_slices(image, buffer=0):
    if len(image.shape) != 3:
        raise ValueError("Image must be a 3D stack")
    focus = stack.compute_focus(image, neighborhood_size=31).mean(axis=(1, 2))
    most_infocus_slice = focus.argmax()
    most_infocus_slices = [most_infocus_slice]
    for i in range(buffer):
        if most_infocus_slice - (i + 1) >= 0:
            most_infocus_slices.append(most_infocus_slice - (i + 1))
        if most_infocus_slice + (i + 1) < image.shape[0]:
            most_infocus_slices.append(most_infocus_slice + (i + 1))
    most_infocus_slices.sort()
    return most_infocus_slices

def main(args):
    # Setup logging
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"extract_most_infocus_{args.Date}_{args.Transcript}_{args.Condition}_{current_time}.log"
    log_dir = os.path.join(args.src_dir, "logs", "extract_most_infocus")
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
    pattern = re.compile(args.Date + "~" + args.Transcript + "~" + args.Condition + "~(.*)" + args.focus_chan + "\.tif")
    for file_name in os.listdir(stack_dir):
        if pattern.search(file_name):
            fov = file_name.replace(f'~{args.focus_chan}.tif', '')
            logging.info(f"Filtering {fov}")
            file = os.path.join(stack_dir, file_name)
            im_stack = stack.read_image(file)
            infocus_slices = pick_most_infocus_slices(im_stack, args.focus_buffer)
            logging.info(f"Most in-focus slices: {infocus_slices}")
            for chan in args.channel_list.split(","):
                file = os.path.join(stack_dir, fov + f'~{chan}.tif')
                im_stack = stack.read_image(file)
                img_filt = im_stack[infocus_slices]
                out_file = file.replace('.tif', f'_filtfocus.tif')
                tifffile.imwrite(out_file, img_filt)
                logging.info(f"Saved filtered image to {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim stacks to just the most in-focus slices.')
    parser.add_argument('Date', type=str, help='Date of the experiment.')
    parser.add_argument('Transcript', type=str, help='Transcript name.')
    parser.add_argument('Condition', type=str, help='Experimental condition.')
    parser.add_argument('src_dir', type=str, help='Source directory containing the data.')
    parser.add_argument('channel_list', type=str, help='List of channels, comma-separated.')
    parser.add_argument('focus_chan', type=str, help='Channel used for focusing.')
    parser.add_argument('focus_buffer', type=int, help='Number of slices to consider around the most in-focus slice.')
    args = parser.parse_args()

    main(args)
