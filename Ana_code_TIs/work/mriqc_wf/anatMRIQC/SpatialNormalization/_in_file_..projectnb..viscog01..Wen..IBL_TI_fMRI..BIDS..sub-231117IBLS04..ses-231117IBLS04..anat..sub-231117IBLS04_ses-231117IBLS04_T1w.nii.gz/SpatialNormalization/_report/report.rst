Node: anatMRIQC (SpatialNormalization (SpatialNormalization (reportlets)
========================================================================


 Hierarchy : mriqc_wf.anatMRIQC.SpatialNormalization.SpatialNormalization
 Exec ID : SpatialNormalization.a0


Original Inputs
---------------


* compress_report : auto
* explicit_masking : True
* flavor : testing
* float : False
* initial_moving_transform : <undefined>
* lesion_mask : <undefined>
* moving : T1w
* moving_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/post_n4/clipped_corrected.nii.gz
* moving_mask : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* num_threads : 1
* orientation : RAS
* out_report : report.svg
* reference : T1w
* reference_image : <undefined>
* reference_mask : /usr2/postdoc/wenwen/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz
* settings : <undefined>
* template : MNI152NLin2009cAsym
* template_resolution : <undefined>
* template_spec : <undefined>


Execution Inputs
----------------


* compress_report : auto
* explicit_masking : True
* flavor : testing
* float : False
* initial_moving_transform : <undefined>
* lesion_mask : <undefined>
* moving : T1w
* moving_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/post_n4/clipped_corrected.nii.gz
* moving_mask : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* num_threads : 1
* orientation : RAS
* out_report : report.svg
* reference : T1w
* reference_image : <undefined>
* reference_mask : /usr2/postdoc/wenwen/.cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz
* settings : <undefined>
* template : MNI152NLin2009cAsym
* template_resolution : <undefined>
* template_spec : <undefined>


Execution Outputs
-----------------


* composite_transform : <undefined>
* elapsed_time : <undefined>
* forward_invert_flags : <undefined>
* forward_transforms : <undefined>
* inverse_composite_transform : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/SpatialNormalization/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/SpatialNormalization/ants_t1_to_mniInverseComposite.h5
* inverse_warped_image : <undefined>
* metric_value : <undefined>
* out_report : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/SpatialNormalization/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/SpatialNormalization/report.svg
* reference_image : <undefined>
* reverse_forward_invert_flags : <undefined>
* reverse_forward_transforms : <undefined>
* reverse_invert_flags : <undefined>
* reverse_transforms : <undefined>
* save_state : <undefined>
* warped_image : <undefined>


Runtime info
------------


* duration : 119.054691
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/SpatialNormalization/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/SpatialNormalization


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
* KMP_DUPLICATE_LIB_OK : True
* KMP_INIT_AT_FORK : FALSE
* LANG : en_US.UTF-8
* LC_ALL : en_US.UTF-8
* LD_LIBRARY_PATH : /usr/lib/x86_64-linux-gnu:/opt/conda/lib:/opt/fsl:/.singularity.d/libs
* MKL_NUM_THREADS : 1
* NIPYPE_NO_ET : 1
* NO_ET : 1
* NSLOTS : 16
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

