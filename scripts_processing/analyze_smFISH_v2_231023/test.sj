#!/bin/bash
set -x  # Enable debugging output

# ... (other parts of the script)

# Loop through each sub-folder in the main directory.
for SUBFOLDER in "${MASK_DIR}"/*; do
    if [ -d "${SUBFOLDER}" ]; then
        # Run the combine_masks.py script for the sub-folder.
        mamba run -n ${ANALYSIS_ENV} python ${combine_masks_prog} "${SUBFOLDER}"

        # Name of the combined mask to copy.
        MASK_NAME=$(basename "${SUBFOLDER}")"_masks.tiff"

        # Copy the combined mask to the combined_masks directory.
        cp "${SUBFOLDER}/${MASK_NAME}" "${COMBINED_MASKS_DIR}"
    fi
done