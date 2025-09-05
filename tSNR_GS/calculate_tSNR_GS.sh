#!/bin/bash
################################################################ 
# This script will calculate the tSNR and GS for OLM participants
data_folder="/home1/Jaquent/Datasets/VP00157_HCP/sessions/"

# .txt listing all the paths
GS_files="GS_paths.txt"
tSNR_files="tSNR_paths.txt"

################################################################ 
# Find all corresponding *dt.series.nii images to be used in the calculation
mapfile -t found_images < <(find ${data_folder} \( -type f -o -type l \) -name "*_AP_Atlas_MSMAll_hp0_clean.dtseries.nii" 2>/dev/null)

# Print to console
echo "##################################################"
echo "${#found_images[@]} images were found." 

################################################################ 
# Loop through the files
for input_file in "${found_images[@]}"; do
    echo "Processing $input_file"

    # Calculate GS and add path to list
    output_file="${input_file/.dtseries.nii/_GScalc.dscalar.nii}"
    wb_command -cifti-reduce $input_file MEAN $output_file
    echo "${output_file}" >> $GS_files

    # Calculate tSNR and add path to list
    output_file="${input_file/.dtseries.nii/_TSNRcalc.dscalar.nii}"
    wb_command -cifti-reduce $input_file TSNR $output_file
    echo "${output_file}" >> $tSNR_files
done

echo "Processing completed"