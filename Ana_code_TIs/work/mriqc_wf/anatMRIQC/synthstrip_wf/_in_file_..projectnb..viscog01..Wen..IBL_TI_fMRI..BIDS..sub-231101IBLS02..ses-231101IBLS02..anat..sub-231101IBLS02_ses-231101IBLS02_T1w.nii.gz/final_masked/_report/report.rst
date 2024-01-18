Node: anatMRIQC (synthstrip_wf (final_masked (nibabel)
======================================================


 Hierarchy : mriqc_wf.anatMRIQC.synthstrip_wf.final_masked
 Exec ID : final_masked.a0


Original Inputs
---------------


* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/post_n4/clipped_corrected.nii.gz
* in_mask : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* threshold : 0.5


Execution Inputs
----------------


* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/post_n4/clipped_corrected.nii.gz
* in_mask : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* threshold : 0.5


Execution Outputs
-----------------


* out_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/final_masked/clipped_corrected_masked.nii.gz


Runtime info
------------


* duration : 0.635823
* hostname : scc-xi1
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/final_masked


Environment
~~~~~~~~~~~


* AFNI_DIR : /opt/afni
* AFNI_IMSAVE_WARNINGS : NO
* AFNI_MODELPATH : /opt/afni/models
* AFNI_PLUGINPATH : /opt/afni/plugins
* AFNI_TTATLAS_DATASET : /opt/afni/atlases
* ANTSPATH : /opt/ants
* CONDA_PATH : /opt/conda
* CPATH : /opt/conda/include:
* FREESURFER_HOME : /opt/freesurfer
* FSLDIR : /opt/fsl
* FSLGECUDAQ : cuda.q
* FSLLOCKDIR : 
* FSLMACHINELIST : 
* FSLMULTIFILEQUIT : TRUE
* FSLOUTPUTTYPE : NIFTI_GZ
* FSLREMOTECALL : 
* FSLTCLSH : /opt/fsl/bin/fsltclsh
* FSLWISH : /opt/fsl/bin/fslwish
* HOME : /usr2/postdoc/wenwen
* IS_DOCKER_8395080871 : 1
* LANG : en_US.UTF-8
* LC_ALL : en_US.UTF-8
* LD_LIBRARY_PATH : /usr/lib/x86_64-linux-gnu:/opt/conda/lib:/opt/fsl:/.singularity.d/libs
* MKL_NUM_THREADS : 1
* NIPYPE_NO_ET : 1
* NO_ET : 1
* NSLOTS : 8
* OMP_NUM_THREADS : 1
* PATH : /opt/fsl/bin:/opt/ants:/opt/afni:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
* POSSUMDIR : /opt/fsl
* PROMPT_COMMAND : PS1="Singularity> "; unset PROMPT_COMMAND
* PS1 : Singularity> 
* PWD : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* PYTHONNOUSERSITE : 1
* PYTHONWARNINGS : ignore
* SINGULARITY_BIND : /share,/usr1,/usr2,/usr3,/usr4,/var/spool/sge,/project,/projectnb,/projectnb2,/restricted,/rproject,/rprojectnb,/rprojectnb2,/scratch,/net,/ad,/var/lib/dbus/machine-id
* SINGULARITY_COMMAND : run
* SINGULARITY_CONTAINER : /share/pkg.7/mriqc/22.0.6/install/bin/mriqc_22.0.6.simg
* SINGULARITY_ENVIRONMENT : /.singularity.d/env/91-environment.sh
* SINGULARITY_NAME : mriqc_22.0.6.simg
* TMPDIR : /scratch/2489098.1.onrcc-m256
* USER : wenwen

