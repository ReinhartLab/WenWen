Node: anatMRIQC (ComputeIQMs (metadata (bids)
=============================================


 Hierarchy : mriqc_wf.anatMRIQC.ComputeIQMs.metadata
 Exec ID : metadata.a0


Original Inputs
---------------


* bids_dir : None
* bids_validate : True
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/sub-231101IBLS02/ses-231101IBLS02/anat/sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz


Execution Inputs
----------------


* bids_dir : None
* bids_validate : True
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/sub-231101IBLS02/ses-231101IBLS02/anat/sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz


Execution Outputs
-----------------


* acquisition : <undefined>
* out_dict : {'AcquisitionMatrixPE': 240, 'AcquisitionNumber': 1, 'AcquisitionTime': '11:55:56.062500', 'BaseResolution': 256, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231101_IBL_S02', 'subject_id': 'CNC_Archive_S02521', 'experiment': '231101_IBL_S02', 'experiment_id': 'CNC_Archive_E02734', 'scan': 7}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3e-06, 'EchoTime': 0.00169, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.00886408, 0.997264, -0.0733949, 0.021191, -0.0735686, -0.997065], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 3, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 93.75, 'PercentSampling': 100, 'PhaseEncodingSteps': 238, 'PhaseResolution': 1, 'PixelBandwidth': 650, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p3', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7;NC1', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 240, 'RefLinesPE': 32, 'RepetitionTime': 2.53, 'SAR': 0.0370875, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_16ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p3 RMS', 'SeriesNumber': 7, 'ShimSetting': [-925, -7137, 3764, 139, -99, -646, -157, -21], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 240.054, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot'}
* reconstruction : <undefined>
* run : <undefined>
* session : 231101IBLS02
* subject : 231101IBLS02
* suffix : <undefined>
* task : <undefined>


Runtime info
------------


* duration : 0.65106
* hostname : scc-xk3
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/metadata


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
* TMPDIR : /scratch/2657761.1.onrcc-m256
* USER : wenwen

