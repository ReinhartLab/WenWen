#!/bin/bash -l

module load fmriprep/22.0.2

# fmriprep command

# bids_dir
# The root folder of a BIDS valid dataset (sub-XXXXX folders should be found at the top level in this folder).
# output_dir
# The output path for the outcomes of preprocessing and visual reports


scratch_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/Prepro/scratch"
bids_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS"
fs_subs_dir="/projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/freesurfer"

# Read the content of the file into a variable
subj=$(cat "$bids_dir/subj_single.txt")
# subj="S02"
sub_label=$(find $bids_dir -maxdepth 1 -type d -name "*sub*$subj")
sub_label=$(basename "$sub_label")

omp=4
mem=16
nproc=4
ncpu=8
dummyN=2 # number of dummy scans, it's the same for localizer task and main task
fd_thresh=0.9 # in mm, following mumford youtube video
dvar_thresh=1.5 #defaults are FD > 0.5 mm or DVARS > 1.5

fmriprep $bids_dir $bids_dir/derivatives/fmriprep participant --participant-label $sub_label -w $scratch_dir/$subj --nprocs $nproc --n-cpus $ncpu --omp-nthreads $omp --mem $mem --use-aroma --output-spaces T1w MNI152NLin2009cAsym:res-2 fsaverage --fd-spike-threshold $fd_thresh --dummy-scans $dummyN --cifti-output 91k --no-submm-recon --dvars-spike-threshold $dvar_thresh --stop-on-first-crash --use-syn-sdc warn --fs-license-file /share/pkg.7/freesurfer/7.3.2/install/license --fs-subjects-dir $fs_subs_dir --skip-bids-validation 


# fmriprep notes:
# https://fmriprep.org/en/20.2.3/outputs.html

# --output-spaces:Standard and non-standard spaces to resample anatomical and functional images to. Standard spaces may be specified by the form <SPACE>[:cohort-<label>][:res-<resolution>]

# For a more accurate estimation of head-motion, we calculate its parameters before any time-domain filtering (i.e., slice-timing correction), as recommended in [Power2017]

# the SliceTiming field is available within the input metadata, this workflow performs slice time correction prior to other signal resampling processes. Slice time correction is performed using AFNI 3dTShift. All slices are realigned in time to the middle of each TR.

# To use ICA-AROMA, set MNI152NLin6Asym:res-2 as a target output space. MELODIC and ICA-AROMA can be run on the resulting images in a separate pipeline.

#fMRIPrep does high-pass filtering before running anatomical or temporal CompCor. Therefore, when using CompCor regressors, the corresponding cosine_XX regressors should also be included in the design matrix.

# If you are already using AROMA-cleaned data (~desc-smoothAROMAnonaggr_bold.nii.gz), do not include ICA-AROMA confounds during your design specification or denoising procedure.

# Nonsteady-states (or dummy scans) in the beginning of every run are dropped before ICA-AROMA is performed. Therefore, any subsequent analysis of ICA-AROMA outputs must drop the same number of nonsteady-states.

# The --dummy-scans sets the number of volumes, but we never remove volumes. The non_steady_state_outlier confound can be used in regression models to keep the dummy volumes from contributing to parameter estimates.

