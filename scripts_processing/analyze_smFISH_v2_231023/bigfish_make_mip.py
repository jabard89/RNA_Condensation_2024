import argparse
import os
import re
import tifffile
import bigfish.stack as stack
import logging
import sys
from datetime import datetime
import matplotlib
import tifffile as tif
matplotlib.use('Agg')

def main(args):
    # Setup logging
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"bigfish_make_mip_{args.Date}_{args.Transcript}_{args.Condition}_{current_time}.log"
    log_dir = os.path.join(args.src_dir, "logs", "bigfish_make_mip")
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
    for trans in os.listdir(date_dir):
        if trans != args.Transcript:
            continue
        trans_dir = os.path.join(date_dir, trans)
        for cond in os.listdir(trans_dir):
            if cond != args.Condition:
                continue
            cond_dir = os.path.join(trans_dir, cond)
            stack_dir = os.path.join(cond_dir, "stacks")
            mip_dir = os.path.join(cond_dir, "mips")
            os.makedirs(mip_dir, exist_ok=True)
            for chan in args.channel_list.split(","):
                for file_name in os.listdir(stack_dir):
                    if args.filt_focus:
                        pattern = re.compile(args.Date + "~" + args.Transcript + "~" + args.Condition + "~(.*)" + chan + "_filtfocus\.tif")
                    else:
                        pattern = re.compile(args.Date + "~" + args.Transcript + "~" + args.Condition + "~(.*)" + chan + "\.tif")
                    if pattern.search(file_name):
                        file = os.path.join(stack_dir, file_name)
                        im_stack = tif.imread(file)
                        im_mip = stack.maximum_projection(im_stack)
                        mip_file = os.path.join(mip_dir, file_name.replace(".tif", "_mip.tif"))
                        tifffile.imwrite(mip_file, im_mip)
                        logging.info(f"Wrote {mip_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create maximum-intensity projections.')
    parser.add_argument('Date', type=str, help='Date of the experiment.')
    parser.add_argument('Transcript', type=str, help='Transcript name.')
    parser.add_argument('Condition', type=str, help='Experimental condition.')
    parser.add_argument('src_dir', type=str, help='Source directory containing the data.')
    parser.add_argument('channel_list', type=str, help='List of channels, comma-separated.')
    parser.add_argument('--filt_focus', action='store_true', default=False, help='Whether to use filtered images.')
    args = parser.parse_args()

    main(args)
