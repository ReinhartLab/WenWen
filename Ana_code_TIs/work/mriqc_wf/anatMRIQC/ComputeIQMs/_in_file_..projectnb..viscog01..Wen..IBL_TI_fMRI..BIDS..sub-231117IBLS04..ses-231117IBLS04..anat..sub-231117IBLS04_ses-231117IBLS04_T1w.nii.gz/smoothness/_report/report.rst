Node: anatMRIQC (ComputeIQMs (smoothness (afni)
===============================================


 Hierarchy : mriqc_wf.anatMRIQC.ComputeIQMs.smoothness
 Exec ID : smoothness.a0


Original Inputs
---------------


* acf : False
* args : -ShowMeClassicFWHM
* arith : <undefined>
* automask : False
* combine : True
* compat : <undefined>
* demed : <undefined>
* detrend : True
* environ : {}
* geom : <undefined>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz
* mask : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* out_detrend : <undefined>
* out_file : <undefined>
* out_subbricks : <undefined>
* unif : <undefined>


Execution Inputs
----------------


* acf : False
* args : -ShowMeClassicFWHM
* arith : <undefined>
* automask : False
* combine : True
* compat : <undefined>
* demed : <undefined>
* detrend : True
* environ : {}
* geom : <undefined>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz
* mask : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz
* out_detrend : <undefined>
* out_file : <undefined>
* out_subbricks : <undefined>
* unif : <undefined>


Execution Outputs
-----------------


* acf_param : <undefined>
* fwhm : (4.21469, 4.06536, 3.90834, 4.06087)
* out_acf : <undefined>
* out_detrend : <undefined>
* out_file : <undefined>
* out_subbricks : <undefined>


Runtime info
------------


* cmdline : 3dFWHMx -ShowMeClassicFWHM -combine -detrend -input /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz -mask /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz -detprefix sub-231117IBLS04_ses-231117IBLS04_T1w_conformed_detrend -out sub-231117IBLS04_ses-231117IBLS04_T1w_conformed_subbricks.out > sub-231117IBLS04_ses-231117IBLS04_T1w_conformed_fwhmx.out
* duration : 6.252821
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/smoothness


Terminal output
~~~~~~~~~~~~~~~


 


Terminal - standard output
~~~~~~~~~~~~~~~~~~~~~~~~~~


 


Terminal - standard error
~~~~~~~~~~~~~~~~~~~~~~~~~


 ++ 3dFWHMx: AFNI version=AFNI_22.2.07 (Aug 19 2022) [64-bit]
++ Authored by: The Bob
[7m*+ WARNING:[0m Using the 'Classic' Gaussian FWHM is not recommended :(
 +  The '-acf' method gives a FWHM estimate which is more robust;
 +  however, assuming the spatial correlation of FMRI noise has
 +  a Gaussian shape is not a good model.
[7m*+ WARNING:[0m   If you are performing spatial transformations on an oblique dset,
  such as /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz,
  or viewing/combining it with volumes of differing obliquity,
  you should consider running: 
     3dWarp -deoblique 
  on this and  other oblique datasets in the same session.
 See 3dWarp -help for details.
++ Oblique dataset:/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz is 5.051599 degrees from plumb.
++ Oblique dataset:/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz is 5.051599 degrees from plumb.
++ Number of voxels in mask = 1390870
[7m*+ WARNING:[0m -demed and/or -corder and/or -unif ignored: only 1 input sub-bricks
++ start Classic FWHM calculations
 + Classic FWHM done (0.00 CPU s thus far)
++ start ACF calculations out to radius = 12.18 mm
 + ACF done (0.00 CPU s thus far)
++ ACF 1D file [radius ACF mixed_model gaussian_NEWmodel] written to 3dFWHMx.1D
++ 1dplot: AFNI version=AFNI_22.2.07 (Aug 19 2022) [64-bit]
++ Authored by: RWC et al.
pnmtopng: 39 colors found
 + and 1dplot-ed to file 3dFWHMx.1D.png


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

