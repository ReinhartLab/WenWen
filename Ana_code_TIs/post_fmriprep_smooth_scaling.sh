#!/bin/bash

module load fsl

bids_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS"
subj=$(cat "$bids_dir/subj_single.txt")
sub_label=$(find $bids_dir -maxdepth 1 -type d -name "*sub*$subj")
sub_label=$(basename "$sub_label")
ses_label=$(find $bids_dir/$sub_label -maxdepth 1 -type d -name "*ses*$subj")
ses_label=$(basename "$ses_label")

func_dir="$bids_dir/derivatives/$sub_label/$ses_label/func"
echo $func_dir

ss_FWHM=4 #spatial smoothing parameter, gaussian kernel in mm, size of 2 voxels
sigma=$(awk "BEGIN { print $ss_FWHM / (sqrt(8 * log(2))) }")

sufix_Str="*space-T1w_desc-preproc_bold.nii.gz" # using T1w space

find "${func_dir}" -name $sufix_Str -type f | while read -r file; do
    # Get the filename without extension
    base_name=$(basename "${file}" .nii.gz)

    echo "$subj: Spatial smoothing with FWHM ${ss_FWHM} + scaling"
    echo $base_name

    # Set the output file path
    output_file="${func_dir}/${base_name}_smooth_scaled.nii.gz"

    # Perform spatial smoothing using fslmaths
    fslmaths "${file}" -s "${sigma}" "${output_file}"

    # Perform scaling using flirt to resample to the desired output resolution
    flirt -in $output_file -ref $output_file -applyisoxfm $output_resolution -out $output_file

    echo "Done. Output: ${output_file}"
done