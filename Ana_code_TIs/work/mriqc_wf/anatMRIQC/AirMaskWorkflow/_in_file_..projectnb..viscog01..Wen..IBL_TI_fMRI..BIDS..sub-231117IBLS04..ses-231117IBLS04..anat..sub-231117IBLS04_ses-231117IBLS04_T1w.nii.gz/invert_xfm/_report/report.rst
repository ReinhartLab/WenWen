Node: anatMRIQC (AirMaskWorkflow (invert_xfm (ants)
===================================================


 Hierarchy : mriqc_wf.anatMRIQC.AirMaskWorkflow.invert_xfm
 Exec ID : invert_xfm.a0


Original Inputs
---------------


* args : <undefined>
* default_value : 0.0
* dimension : 3
* environ : {'NSLOTS': '1'}
* float : True
* input_image : /usr2/postdoc/wenwen/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_desc-head_mask.nii.gz
* input_image_type : <undefined>
* interpolation : MultiLabel
* interpolation_parameters : <undefined>
* invert_transform_flags : <undefined>
* num_threads : 1
* out_postfix : _trans
* output_image : <undefined>
* print_out_composite_warp_file : <undefined>
* reference_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* transforms : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/SpatialNormalization/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/SpatialNormalization/ants_t1_to_mniInverseComposite.h5']


Execution Inputs
----------------


* args : <undefined>
* default_value : 0.0
* dimension : 3
* environ : {'NSLOTS': '1'}
* float : True
* input_image : /usr2/postdoc/wenwen/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_desc-head_mask.nii.gz
* input_image_type : <undefined>
* interpolation : MultiLabel
* interpolation_parameters : <undefined>
* invert_transform_flags : <undefined>
* num_threads : 1
* out_postfix : _trans
* output_image : <undefined>
* print_out_composite_warp_file : <undefined>
* reference_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* transforms : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/SpatialNormalization/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/SpatialNormalization/ants_t1_to_mniInverseComposite.h5']


Execution Outputs
-----------------


* output_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/invert_xfm/tpl-MNI152NLin2009cAsym_res-01_desc-head_mask_trans.nii.gz


Runtime info
------------


* cmdline : antsApplyTransforms --default-value 0 --dimensionality 3 --float 1 --input /usr2/postdoc/wenwen/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_desc-head_mask.nii.gz --interpolation MultiLabel --output tpl-MNI152NLin2009cAsym_res-01_desc-head_mask_trans.nii.gz --reference-image /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz --transform /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/SpatialNormalization/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/SpatialNormalization/ants_t1_to_mniInverseComposite.h5
* duration : 85.898677
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/invert_xfm


Terminal output
~~~~~~~~~~~~~~~


 


Terminal - standard output
~~~~~~~~~~~~~~~~~~~~~~~~~~


 


Terminal - standard error
~~~~~~~~~~~~~~~~~~~~~~~~~


 


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
* NSLOTS : 1
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
* TMPDIR : /scratch/2657782.1.onrcc-m256
* USER : wenwen

