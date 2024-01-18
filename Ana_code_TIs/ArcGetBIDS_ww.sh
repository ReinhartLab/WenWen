#!/bin/bash

read -p "Enter ScanName (e.g., 231016_IBL_S01): " ScanName

read -p "Use customized yaml for configuration (y/n)? " response

if [ "${response,,}" == "n" ]; then
  yaml="IBL_general"
else
  read -p "Which yaml file for this subject [e.g., IBL_S01] ): " yaml
fi

download_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS"
project_name="Reinhart_IBL"
yaml_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/yamls"
yaml_file="$yaml_dir/$yaml".yaml""

# Print the download directory and project name

echo "Download directory: $download_dir"
echo "Project: $project_name"
echo "yaml file: $yaml_file"

# Load modules
module load python3
module load yaxil/0.6.5

# Get data
ArcGet.py -a xnat -s "$ScanName" -p "$project_name" -f bids -c $yaml_file -o "$download_dir"

# sourcedata of each participant will be automatically moved to a shared folder