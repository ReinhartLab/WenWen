
#!/bin/bash -l

module load mriqc/22.0.6

# reference link at https://mriqc.readthedocs.io/en/latest/running.html#command-line-interface

bids_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS"

# Read the content of the file into a variable
subj=$(cat "$bids_dir/subj_single.txt")
# subj="S02"

sub_label=$(find $bids_dir -maxdepth 1 -type d -name "*sub*$subj")
sub_label=$(basename "$sub_label")
echo $sub_label

# run mriqc on indicated subjects
mriqc $bids_dir/ $bids_dir/derivatives/mriqc/ participant --participant-label $sub_label

# generate summary for mriqc output under specified directory
mriqc $bids_dir/ $bids_dir/derivatives/mriqc/ group

# group aggregation:
# each bold scan is one dot
# click matric name for explanations
