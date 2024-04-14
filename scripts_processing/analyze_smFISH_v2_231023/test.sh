#!/bin/bash
# Jared Bard
# 231020
# pipeline for analyzing smFISH data
# 1. extract most in focus slices (here using top 3)
# 2. make maximum intensity projections (mips) for those three slices
# 3. go to cell profiler, and use make_mip_masks.cppipe to extract masks
# 4. combine masks together
# 5. run bigfish on the mips
# 6. analyze the results
# 7. collate results into a single csv file

# Parameters
SRC_DIR="/Users/Caitlin/Dropbox (Drummond Lab)/smFISH_data/fishquant"
declare -a DATES=("20230810" "20230816" "20230823" "20230825")
declare -a TRANSCRIPTS=("SSA4" "SSB1" "ADD66" "HSP104")
declare -a CONDITIONS=("30C" "42C" "46C")

RNA_CHAN="647"
PAB1_CHAN="A488"
FOCUS_BUFFER="1"
THRESH="400"
VOX_X="130"
VOX_Y="130"
VOX_Z="200"
CHANNEL_LIST="647,A488"
RANDPOINTS="100"
RNA_BUFFER="1"

PROGRAM_DIR="/Users/Caitlin/Dropbox (Drummond Lab)/smFISH_data/analyze_smFISH_v2_231023"
extract_most_infocus_prog="${PROGRAM_DIR}/extract_most_infocus.py"
bigfish_make_mip_prog="${PROGRAM_DIR}/bigfish_make_mip.py"
bigfish_test_thresh_prog="${PROGRAM_DIR}/bigfish_test_thresh.py"
bigfish_run_all_prog="${PROGRAM_DIR}/bigfish_run_all.py"
bigfish_analyze_prog="${PROGRAM_DIR}/analyze_smFISH.py"
combine_masks_prog="${PROGRAM_DIR}/combine_masks.py"
collate_results_script="${PROGRAM_DIR}/collate_smFISH_analysis_outputs.R"

BIGFISH_ENV="bigfish_env"
CELLPOSE_ENV="cellpose"
ANALYSIS_ENV="images"

# # Extract most in-focus and make MIPs
# for DATE in "${DATES[@]}"; do
#     for TRANSCRIPT in "${TRANSCRIPTS[@]}"; do
#         for CONDITION in "${CONDITIONS[@]}"; do
#             cond_dir="${SRC_DIR}/${DATE}/${TRANSCRIPT}/${CONDITION}"
#             if [ -d "${cond_dir}/stacks" ]; then
#                 echo "${cond_dir}"
#                 conda run -n ${BIGFISH_ENV} python "${extract_most_infocus_prog}" "${DATE}" "${TRANSCRIPT}" "${CONDITION}" "${SRC_DIR}" "${CHANNEL_LIST}" "${PAB1_CHAN}" "${FOCUS_BUFFER}"
#                 conda run -n ${BIGFISH_ENV} python "${bigfish_make_mip_prog}" "${DATE}" "${TRANSCRIPT}" "${CONDITION}" "${SRC_DIR}" "647,A488" --filt_focus
#             fi
#         done
#     done
# done

# # Run CellProfiler make_mip_masks.cppipe here

# # Combine masks
# MASK_DIR="${SRC_DIR}/masks"
# COMBINED_MASKS_DIR="${SRC_DIR}/combined_masks"

# # Create the combined_masks directory if it doesn't exist.
# mkdir -p "${COMBINED_MASKS_DIR}"

# # Loop through each sub-folder in the main directory.
# for SUBFOLDER in "${MASK_DIR}"/*; do
#     if [ -d "${SUBFOLDER}" ]; then
#         conda run -n ${ANALYSIS_ENV} python "${combine_masks_prog}" "${SUBFOLDER}"

#         # Name of the combined mask to copy.
#         MASK_NAME=$(basename "${SUBFOLDER}")"_masks.tiff"

#         # Copy the combined mask to the combined_masks directory.
#         cp "${SUBFOLDER}/${MASK_NAME}" "${COMBINED_MASKS_DIR}"
#     fi
# done

# # Run bigfish_run_all.py, analyze_smFISH.py, and collate results
# for DATE in "${DATES[@]}"; do
#     for TRANSCRIPT in "${TRANSCRIPTS[@]}"; do
#         for CONDITION in "${CONDITIONS[@]}"; do
#             cond_dir="${SRC_DIR}/${DATE}/${TRANSCRIPT}/${CONDITION}"
#             if [ -d "${cond_dir}/stacks" ]; then
#                 echo "${cond_dir}"
#                 conda run -n ${BIGFISH_ENV} python "${bigfish_run_all_prog}" "${DATE}" "${TRANSCRIPT}" "${CONDITION}" "${SRC_DIR}" "${RNA_CHAN}" "${VOX_X}" "${VOX_Y}" "${VOX_Z}" "${THRESH}" --filt_focus
#                 conda run -n ${ANALYSIS_ENV} python "${bigfish_analyze_prog}" "${DATE}" "${TRANSCRIPT}" "${CONDITION}" "${SRC_DIR}" "${RNA_CHAN}" "${PAB1_CHAN}" "${THRESH}" "${RANDPOINTS}" "${RNA_BUFFER}" --filt_focus
#             fi
#         done
#     done
# done

# Collate results
Rscript --vanilla "${collate_results_script}"

echo "All programs completed!"
