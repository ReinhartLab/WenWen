#!/bin/bash

# Prompt for ScanName input
read -p "Enter ScanName (e.g., 231016_IBL_S01): " ScanName

# Prompt for SubjID input
read -p "Enter SubjID (e.g., S01): " SubjID

# Directory where the data will be downloaded
download_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/Raw/$SubjID"
project_name="Reinhart_IBL"

# Print the download directory and project name

echo "Download directory: $download_dir"
echo "Project: $project_name"

# Load modules
module load python3
module load yaxil/0.6.5

# Get data
ArcGet.py -a xnat -s "$ScanName" -p "$project_name" -f native -o "$download_dir"

# move all scans from the "scans" folder and delete that folder
export subj_dir_orig=${download_dir}/$ScanName/
export subj_dir_new=${download_dir}/

echo "$subj_dir_new"
echo "$subj_dir_orig"

mv $subj_dir_orig/scans/* $subj_dir_new
rm -r $subj_dir_orig

# loop through scans folders and remove annoying nested directories (adapted from Jasmine_Ling_Lab)
 for dir in ${subj_dir_new}/*
 do
 	echo "$dir"
 	mv {$dir}/resources/DICOM/files/* $dir
 	rm -r $dir/resources/
 done

 # interactive prompt: 
read -p "Do you need to delete some scans (Y/N)? " response

if [ "${response^^}" == "Y" ]; then
  read -p "Enter scan ID separated by comma (e.g., 27,28): " input

  # Split the comma-separated input into an array
  IFS=',' read -ra names <<< "$input"

  for folder in "$subj_dir_new"/*; do
    if [ -d "$folder" ]; then
      folder_name=$(basename "$folder")

      for name in "${names[@]}"; do

        case "$folder_name" in
          *"$name"-""*)
          # Ask for confirmation before deleting
          read -p "Delete folder '$folder_name'? (Y/N): " confirm
          if [ "${confirm^^}" == "Y" ]; then
              # Delete the folder
              rm -r "$folder"
              echo "Folder '$folder_name' deleted."
            else
              echo "Folder '$folder_name' not deleted."
            fi
            ;;
        esac

      done
    fi
	done
fi

echo "Completed!!!"

