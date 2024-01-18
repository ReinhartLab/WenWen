Node: anatMRIQC (HeadMaskWorkflow (fsl_bet (fsl)
================================================


 Hierarchy : mriqc_wf.anatMRIQC.HeadMaskWorkflow.fsl_bet
 Exec ID : fsl_bet.a0


Original Inputs
---------------


* args : <undefined>
* center : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* frac : <undefined>
* functional : <undefined>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/HeadMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/fsl_bet/clipped_corrected.nii.gz
* mask : <undefined>
* mesh : <undefined>
* no_output : <undefined>
* out_file : <undefined>
* outline : <undefined>
* output_type : NIFTI_GZ
* padding : <undefined>
* radius : <undefined>
* reduce_bias : <undefined>
* remove_eyes : <undefined>
* robust : <undefined>
* skull : <undefined>
* surfaces : True
* t2_guided : <undefined>
* threshold : <undefined>
* vertical_gradient : <undefined>


Execution Inputs
----------------


* args : <undefined>
* center : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* frac : <undefined>
* functional : <undefined>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/HeadMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/fsl_bet/clipped_corrected.nii.gz
* mask : <undefined>
* mesh : <undefined>
* no_output : <undefined>
* out_file : <undefined>
* outline : <undefined>
* output_type : NIFTI_GZ
* padding : <undefined>
* radius : <undefined>
* reduce_bias : <undefined>
* remove_eyes : <undefined>
* robust : <undefined>
* skull : <undefined>
* surfaces : True
* t2_guided : <undefined>
* threshold : <undefined>
* vertical_gradient : <undefined>


Execution Outputs
-----------------


* inskull_mask_file : <undefined>
* inskull_mesh_file : <undefined>
* mask_file : <undefined>
* meshfile : <undefined>
* out_file : <undefined>
* outline_file : <undefined>
* outskin_mask_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/HeadMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/fsl_bet/clipped_corrected_brain_outskin_mask.nii.gz
* outskin_mesh_file : <undefined>
* outskull_mask_file : <undefined>
* outskull_mesh_file : <undefined>
* skull_file : <undefined>
* skull_mask_file : <undefined>


Runtime info
------------


* cmdline : bet clipped_corrected.nii.gz clipped_corrected_brain.nii.gz -A
* duration : 59.409167
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/HeadMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/fsl_bet


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

