Node: anatMRIQC (segmentation (fsl)
===================================


 Hierarchy : mriqc_wf.anatMRIQC.segmentation
 Exec ID : segmentation.a0


Original Inputs
---------------


* args : <undefined>
* bias_iters : <undefined>
* bias_lowpass : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* hyper : <undefined>
* img_type : 1
* in_files : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/clipped_corrected_masked.nii.gz']
* init_seg_smooth : <undefined>
* init_transform : <undefined>
* iters_afterbias : <undefined>
* manual_seg : <undefined>
* mixel_smooth : <undefined>
* no_bias : <undefined>
* no_pve : <undefined>
* number_classes : <undefined>
* other_priors : <undefined>
* out_basename : segment
* output_biascorrected : <undefined>
* output_biasfield : <undefined>
* output_type : NIFTI_GZ
* probability_maps : <undefined>
* segment_iters : <undefined>
* segments : True
* use_priors : <undefined>
* verbose : <undefined>


Execution Inputs
----------------


* args : <undefined>
* bias_iters : <undefined>
* bias_lowpass : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* hyper : <undefined>
* img_type : 1
* in_files : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/clipped_corrected_masked.nii.gz']
* init_seg_smooth : <undefined>
* init_transform : <undefined>
* iters_afterbias : <undefined>
* manual_seg : <undefined>
* mixel_smooth : <undefined>
* no_bias : <undefined>
* no_pve : <undefined>
* number_classes : <undefined>
* other_priors : <undefined>
* out_basename : segment
* output_biascorrected : <undefined>
* output_biasfield : <undefined>
* output_type : NIFTI_GZ
* probability_maps : <undefined>
* segment_iters : <undefined>
* segments : True
* use_priors : <undefined>
* verbose : <undefined>


Execution Outputs
-----------------


* bias_field : <undefined>
* mixeltype : <undefined>
* partial_volume_files : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/segment_pve_0.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/segment_pve_1.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/segment_pve_2.nii.gz']
* partial_volume_map : <undefined>
* probability_maps : <undefined>
* restored_image : <undefined>
* tissue_class_files : <undefined>
* tissue_class_map : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/segment_seg.nii.gz


Runtime info
------------


* cmdline : fast -t 1 -o segment -g -S 1 /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation/clipped_corrected_masked.nii.gz
* duration : 206.510959
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/segmentation


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

