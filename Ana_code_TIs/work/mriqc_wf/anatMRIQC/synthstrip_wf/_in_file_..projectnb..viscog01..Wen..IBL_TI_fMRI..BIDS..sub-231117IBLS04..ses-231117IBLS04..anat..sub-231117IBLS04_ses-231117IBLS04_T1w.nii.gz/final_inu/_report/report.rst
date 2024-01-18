Node: anatMRIQC (synthstrip_wf (final_inu (utility)
===================================================


 Hierarchy : mriqc_wf.anatMRIQC.synthstrip_wf.final_inu
 Exec ID : final_inu.a0


Original Inputs
---------------


* bias_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/post_n4/clipped_bias.nii.gz
* function_str : def _apply_bias_correction(in_file, bias_image, out_file=None):
    import os.path as op

    import numpy as np
    import nibabel as nb

    img = nb.load(in_file)
    data = np.clip(
        img.get_fdata() * nb.load(bias_image).get_fdata(),
        a_min=0,
        a_max=None,
    )
    out_img = img.__class__(
        data.astype(img.get_data_dtype()),
        img.affine,
        img.header,
    )

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_file = op.abspath("{}_inu{}".format(fname, ext))

    out_img.to_filename(out_file)
    return out_file

* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz
* out_file : <undefined>


Execution Inputs
----------------


* bias_image : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/post_n4/clipped_bias.nii.gz
* function_str : def _apply_bias_correction(in_file, bias_image, out_file=None):
    import os.path as op

    import numpy as np
    import nibabel as nb

    img = nb.load(in_file)
    data = np.clip(
        img.get_fdata() * nb.load(bias_image).get_fdata(),
        a_min=0,
        a_max=None,
    )
    out_img = img.__class__(
        data.astype(img.get_data_dtype()),
        img.affine,
        img.header,
    )

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_file = op.abspath("{}_inu{}".format(fname, ext))

    out_img.to_filename(out_file)
    return out_file

* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/conform/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed.nii.gz
* out_file : <undefined>


Execution Outputs
-----------------


* out : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/final_inu/sub-231117IBLS04_ses-231117IBLS04_T1w_conformed_inu.nii.gz


Runtime info
------------


* duration : 0.952369
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/final_inu


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

