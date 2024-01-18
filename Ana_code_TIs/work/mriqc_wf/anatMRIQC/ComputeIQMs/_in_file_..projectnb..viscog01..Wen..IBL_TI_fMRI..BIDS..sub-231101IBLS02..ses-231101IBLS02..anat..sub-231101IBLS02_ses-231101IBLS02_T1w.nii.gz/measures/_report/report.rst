Node: anatMRIQC (ComputeIQMs (measures (anatomical)
===================================================


 Hierarchy : mriqc_wf.anatMRIQC.ComputeIQMs.measures
 Exec ID : measures.a0


Original Inputs
---------------


* air_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/ArtifactMask/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed_air.nii.gz
* artifact_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/ArtifactMask/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed_art.nii.gz
* head_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/HeadMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/fsl_bet/clipped_corrected_brain_outskin_mask.nii.gz
* human : True
* in_bias : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/post_n4/clipped_bias.nii.gz
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/conform/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed.nii.gz
* in_fwhm : [4.26038, 4.14507, 3.76862]
* in_noinu : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/harmonize/clipped_corrected_harmonized.nii.gz
* in_pvms : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_pve_0.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_pve_1.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_pve_2.nii.gz']
* in_segm : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_seg.nii.gz
* in_tpms : <undefined>
* mni_tpms : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/MNItpms2t1/mapflow/_MNItpms2t10/tpl-MNI152NLin2009cAsym_res-01_label-CSF_probseg_trans.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/MNItpms2t1/mapflow/_MNItpms2t11/tpl-MNI152NLin2009cAsym_res-01_label-GM_probseg_trans.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/MNItpms2t1/mapflow/_MNItpms2t12/tpl-MNI152NLin2009cAsym_res-01_label-WM_probseg_trans.nii.gz']
* rot_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/RotationMask/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed_rotmask.nii.gz


Execution Inputs
----------------


* air_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/ArtifactMask/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed_air.nii.gz
* artifact_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/ArtifactMask/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed_art.nii.gz
* head_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/HeadMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/fsl_bet/clipped_corrected_brain_outskin_mask.nii.gz
* human : True
* in_bias : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/synthstrip_wf/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/post_n4/clipped_bias.nii.gz
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/conform/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed.nii.gz
* in_fwhm : [4.26038, 4.14507, 3.76862]
* in_noinu : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/harmonize/clipped_corrected_harmonized.nii.gz
* in_pvms : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_pve_0.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_pve_1.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_pve_2.nii.gz']
* in_segm : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/segmentation/segment_seg.nii.gz
* in_tpms : <undefined>
* mni_tpms : ['/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/MNItpms2t1/mapflow/_MNItpms2t10/tpl-MNI152NLin2009cAsym_res-01_label-CSF_probseg_trans.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/MNItpms2t1/mapflow/_MNItpms2t11/tpl-MNI152NLin2009cAsym_res-01_label-GM_probseg_trans.nii.gz', '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/MNItpms2t1/mapflow/_MNItpms2t12/tpl-MNI152NLin2009cAsym_res-01_label-WM_probseg_trans.nii.gz']
* rot_msk : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/AirMaskWorkflow/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/RotationMask/sub-231101IBLS02_ses-231101IBLS02_T1w_conformed_rotmask.nii.gz


Execution Outputs
-----------------


* cjv : <undefined>
* cnr : <undefined>
* efc : <undefined>
* fber : <undefined>
* fwhm : <undefined>
* icvs : <undefined>
* inu : <undefined>
* out_noisefit : <undefined>
* out_qc : {'summary_csf_mean': 243.7133794456424, 'summary_csf_stdv': 99.32487335290412, 'summary_csf_median': 229.14835105091333, 'summary_csf_mad': 93.72155727603942, 'summary_csf_p95': 435.225447328016, 'summary_csf_p05': 109.49376810342073, 'summary_csf_k': 0.880289115292312, 'summary_csf_n': 55608.0, 'summary_gm_mean': 697.024150145573, 'summary_gm_stdv': 71.26405326432702, 'summary_gm_median': 695.6150722131133, 'summary_gm_mad': 69.41911665722183, 'summary_gm_p95': 818.1188200354576, 'summary_gm_p05': 580.6039557948709, 'summary_gm_k': 0.05697428415867156, 'summary_gm_n': 33149.0, 'summary_wm_mean': 999.3077666205476, 'summary_wm_stdv': 42.36565429495651, 'summary_wm_median': 999.9741367027164, 'summary_wm_mad': 40.44018115591885, 'summary_wm_p95': 1067.2216670066118, 'summary_wm_p05': 928.6157423779368, 'summary_wm_k': 0.41421103827590633, 'summary_wm_n': 237006.0, 'summary_bg_mean': 7.294760387976912, 'summary_bg_stdv': 16.726895529030674, 'summary_bg_median': 4.369408927857876, 'summary_bg_mad': 6.47809537000027, 'summary_bg_p95': 24.97543801367283, 'summary_bg_p05': 0.0, 'summary_bg_k': 1415.698565637623, 'summary_bg_n': 2943587.0, 'snr_csf': 2.3070383370873855, 'snr_wm': 23.603365597562945, 'snr_gm': 9.760946059378423, 'snr_total': 11.890449998009585, 'snrd_csf': 23.174006741031555, 'snrd_wm': 101.1284055875976, 'snrd_gm': 70.34826259364256, 'snrd_total': 64.88355830742391, 'cnr': 3.598619806669339, 'fber': 11288.292585591262, 'efc': 0.5773906520868097, 'wm2max': 0.6809177479722445, 'qi_1': 2.0223060355723633e-05, 'cjv': 0.3609529356300591, 'fwhm_x': 4.26038, 'fwhm_y': 4.14507, 'fwhm_z': 3.76862, 'fwhm_avg': 4.058023333333333, 'icvs_csf': 0.22255566176794614, 'icvs_gm': 0.42319373425092205, 'icvs_wm': 0.35425060398113184, 'rpve_csf': 19.978500611846922, 'rpve_gm': 10.771963312626635, 'rpve_wm': 14.86800622177038, 'size_x': 176, 'size_y': 240, 'size_z': 256, 'spacing_x': 1.0, 'spacing_y': 1.0, 'spacing_z': 1.0, 'inu_range': 0.21127622574567795, 'inu_med': 0.5270328521728516, 'tpm_overlap_csf': 0.16608310909615706, 'tpm_overlap_gm': 0.4676646456480538, 'tpm_overlap_wm': 0.521670114151025}
* qi_1 : <undefined>
* rpve : <undefined>
* size : <undefined>
* snr : <undefined>
* snrd : <undefined>
* spacing : <undefined>
* summary : <undefined>
* tpm_overlap : <undefined>
* wm2max : <undefined>


Runtime info
------------


* duration : 7.510101
* hostname : scc-xj1
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/measures


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
* TMPDIR : /scratch/2489457.1.onrcc-m256
* USER : wenwen

