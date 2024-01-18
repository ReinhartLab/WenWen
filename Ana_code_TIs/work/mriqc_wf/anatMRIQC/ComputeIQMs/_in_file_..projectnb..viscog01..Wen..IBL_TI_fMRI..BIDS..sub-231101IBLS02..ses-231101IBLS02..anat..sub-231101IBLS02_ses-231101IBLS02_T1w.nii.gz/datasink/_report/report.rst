Node: anatMRIQC (ComputeIQMs (datasink (bids)
=============================================


 Hierarchy : mriqc_wf.anatMRIQC.ComputeIQMs.datasink
 Exec ID : datasink.a0


Original Inputs
---------------


* _outputs : {'qi_2': 0.003449920283194974}
* acq_id : <undefined>
* dataset : <unset>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/sub-231101IBLS02/ses-231101IBLS02/anat/sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz
* metadata : {'AcquisitionMatrixPE': 240, 'AcquisitionNumber': 1, 'AcquisitionTime': '11:55:56.062500', 'BaseResolution': 256, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231101_IBL_S02', 'subject_id': 'CNC_Archive_S02521', 'experiment': '231101_IBL_S02', 'experiment_id': 'CNC_Archive_E02734', 'scan': 7}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3e-06, 'EchoTime': 0.00169, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.00886408, 0.997264, -0.0733949, 0.021191, -0.0735686, -0.997065], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 3, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 93.75, 'PercentSampling': 100, 'PhaseEncodingSteps': 238, 'PhaseResolution': 1, 'PixelBandwidth': 650, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p3', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7;NC1', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 240, 'RefLinesPE': 32, 'RepetitionTime': 2.53, 'SAR': 0.0370875, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_16ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p3 RMS', 'SeriesNumber': 7, 'ShimSetting': [-925, -7137, 3764, 139, -99, -646, -157, -21], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 240.054, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot'}
* modality : T1w
* out_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/mriqc
* provenance : {'md5sum': 'c9652695037ed45b6a0df330b0c813c6', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}
* rec_id : <undefined>
* root : {'summary_csf_mean': 243.7133794456424, 'summary_csf_stdv': 99.32487335290412, 'summary_csf_median': 229.14835105091333, 'summary_csf_mad': 93.72155727603942, 'summary_csf_p95': 435.225447328016, 'summary_csf_p05': 109.49376810342073, 'summary_csf_k': 0.880289115292312, 'summary_csf_n': 55608.0, 'summary_gm_mean': 697.024150145573, 'summary_gm_stdv': 71.26405326432702, 'summary_gm_median': 695.6150722131133, 'summary_gm_mad': 69.41911665722183, 'summary_gm_p95': 818.1188200354576, 'summary_gm_p05': 580.6039557948709, 'summary_gm_k': 0.05697428415867156, 'summary_gm_n': 33149.0, 'summary_wm_mean': 999.3077666205476, 'summary_wm_stdv': 42.36565429495651, 'summary_wm_median': 999.9741367027164, 'summary_wm_mad': 40.44018115591885, 'summary_wm_p95': 1067.2216670066118, 'summary_wm_p05': 928.6157423779368, 'summary_wm_k': 0.41421103827590633, 'summary_wm_n': 237006.0, 'summary_bg_mean': 7.294760387976912, 'summary_bg_stdv': 16.726895529030674, 'summary_bg_median': 4.369408927857876, 'summary_bg_mad': 6.47809537000027, 'summary_bg_p95': 24.97543801367283, 'summary_bg_p05': 0.0, 'summary_bg_k': 1415.698565637623, 'summary_bg_n': 2943587.0, 'snr_csf': 2.3070383370873855, 'snr_wm': 23.603365597562945, 'snr_gm': 9.760946059378423, 'snr_total': 11.890449998009585, 'snrd_csf': 23.174006741031555, 'snrd_wm': 101.1284055875976, 'snrd_gm': 70.34826259364256, 'snrd_total': 64.88355830742391, 'cnr': 3.598619806669339, 'fber': 11288.292585591262, 'efc': 0.5773906520868097, 'wm2max': 0.6809177479722445, 'qi_1': 2.0223060355723633e-05, 'cjv': 0.3609529356300591, 'fwhm_x': 4.26038, 'fwhm_y': 4.14507, 'fwhm_z': 3.76862, 'fwhm_avg': 4.058023333333333, 'icvs_csf': 0.22255566176794614, 'icvs_gm': 0.42319373425092205, 'icvs_wm': 0.35425060398113184, 'rpve_csf': 19.978500611846922, 'rpve_gm': 10.771963312626635, 'rpve_wm': 14.86800622177038, 'size_x': 176, 'size_y': 240, 'size_z': 256, 'spacing_x': 1.0, 'spacing_y': 1.0, 'spacing_z': 1.0, 'inu_range': 0.21127622574567795, 'inu_med': 0.5270328521728516, 'tpm_overlap_csf': 0.16608310909615706, 'tpm_overlap_gm': 0.4676646456480538, 'tpm_overlap_wm': 0.521670114151025, 'qi_2': 0.003449920283194974, 'bids_meta': {'subject_id': '231101IBLS02', 'session_id': '231101IBLS02', 'modality': 'T1w', 'AcquisitionMatrixPE': 240, 'AcquisitionNumber': 1, 'AcquisitionTime': '11:55:56.062500', 'BaseResolution': 256, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231101_IBL_S02', 'subject_id': 'CNC_Archive_S02521', 'experiment': '231101_IBL_S02', 'experiment_id': 'CNC_Archive_E02734', 'scan': 7}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3e-06, 'EchoTime': 0.00169, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.00886408, 0.997264, -0.0733949, 0.021191, -0.0735686, -0.997065], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 3, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 93.75, 'PercentSampling': 100, 'PhaseEncodingSteps': 238, 'PhaseResolution': 1, 'PixelBandwidth': 650, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p3', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7;NC1', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 240, 'RefLinesPE': 32, 'RepetitionTime': 2.53, 'SAR': 0.0370875, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_16ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p3 RMS', 'SeriesNumber': 7, 'ShimSetting': [-925, -7137, 3764, 139, -99, -646, -157, -21], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 240.054, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot', 'dataset': '<unset>'}, 'provenance': {'md5sum': 'c9652695037ed45b6a0df330b0c813c6', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}}
* run_id : <undefined>
* session_id : 231101IBLS02
* subject_id : 231101IBLS02
* task_id : <undefined>


Execution Inputs
----------------


* _outputs : {'qi_2': 0.003449920283194974}
* acq_id : <undefined>
* dataset : <unset>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/sub-231101IBLS02/ses-231101IBLS02/anat/sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz
* metadata : {'AcquisitionMatrixPE': 240, 'AcquisitionNumber': 1, 'AcquisitionTime': '11:55:56.062500', 'BaseResolution': 256, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231101_IBL_S02', 'subject_id': 'CNC_Archive_S02521', 'experiment': '231101_IBL_S02', 'experiment_id': 'CNC_Archive_E02734', 'scan': 7}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3e-06, 'EchoTime': 0.00169, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.00886408, 0.997264, -0.0733949, 0.021191, -0.0735686, -0.997065], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 3, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 93.75, 'PercentSampling': 100, 'PhaseEncodingSteps': 238, 'PhaseResolution': 1, 'PixelBandwidth': 650, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p3', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7;NC1', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 240, 'RefLinesPE': 32, 'RepetitionTime': 2.53, 'SAR': 0.0370875, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_16ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p3 RMS', 'SeriesNumber': 7, 'ShimSetting': [-925, -7137, 3764, 139, -99, -646, -157, -21], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 240.054, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot'}
* modality : T1w
* out_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/mriqc
* provenance : {'md5sum': 'c9652695037ed45b6a0df330b0c813c6', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}
* rec_id : <undefined>
* root : {'summary_csf_mean': 243.7133794456424, 'summary_csf_stdv': 99.32487335290412, 'summary_csf_median': 229.14835105091333, 'summary_csf_mad': 93.72155727603942, 'summary_csf_p95': 435.225447328016, 'summary_csf_p05': 109.49376810342073, 'summary_csf_k': 0.880289115292312, 'summary_csf_n': 55608.0, 'summary_gm_mean': 697.024150145573, 'summary_gm_stdv': 71.26405326432702, 'summary_gm_median': 695.6150722131133, 'summary_gm_mad': 69.41911665722183, 'summary_gm_p95': 818.1188200354576, 'summary_gm_p05': 580.6039557948709, 'summary_gm_k': 0.05697428415867156, 'summary_gm_n': 33149.0, 'summary_wm_mean': 999.3077666205476, 'summary_wm_stdv': 42.36565429495651, 'summary_wm_median': 999.9741367027164, 'summary_wm_mad': 40.44018115591885, 'summary_wm_p95': 1067.2216670066118, 'summary_wm_p05': 928.6157423779368, 'summary_wm_k': 0.41421103827590633, 'summary_wm_n': 237006.0, 'summary_bg_mean': 7.294760387976912, 'summary_bg_stdv': 16.726895529030674, 'summary_bg_median': 4.369408927857876, 'summary_bg_mad': 6.47809537000027, 'summary_bg_p95': 24.97543801367283, 'summary_bg_p05': 0.0, 'summary_bg_k': 1415.698565637623, 'summary_bg_n': 2943587.0, 'snr_csf': 2.3070383370873855, 'snr_wm': 23.603365597562945, 'snr_gm': 9.760946059378423, 'snr_total': 11.890449998009585, 'snrd_csf': 23.174006741031555, 'snrd_wm': 101.1284055875976, 'snrd_gm': 70.34826259364256, 'snrd_total': 64.88355830742391, 'cnr': 3.598619806669339, 'fber': 11288.292585591262, 'efc': 0.5773906520868097, 'wm2max': 0.6809177479722445, 'qi_1': 2.0223060355723633e-05, 'cjv': 0.3609529356300591, 'fwhm_x': 4.26038, 'fwhm_y': 4.14507, 'fwhm_z': 3.76862, 'fwhm_avg': 4.058023333333333, 'icvs_csf': 0.22255566176794614, 'icvs_gm': 0.42319373425092205, 'icvs_wm': 0.35425060398113184, 'rpve_csf': 19.978500611846922, 'rpve_gm': 10.771963312626635, 'rpve_wm': 14.86800622177038, 'size_x': 176, 'size_y': 240, 'size_z': 256, 'spacing_x': 1.0, 'spacing_y': 1.0, 'spacing_z': 1.0, 'inu_range': 0.21127622574567795, 'inu_med': 0.5270328521728516, 'tpm_overlap_csf': 0.16608310909615706, 'tpm_overlap_gm': 0.4676646456480538, 'tpm_overlap_wm': 0.521670114151025, 'qi_2': 0.003449920283194974, 'bids_meta': {'subject_id': '231101IBLS02', 'session_id': '231101IBLS02', 'modality': 'T1w', 'AcquisitionMatrixPE': 240, 'AcquisitionNumber': 1, 'AcquisitionTime': '11:55:56.062500', 'BaseResolution': 256, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231101_IBL_S02', 'subject_id': 'CNC_Archive_S02521', 'experiment': '231101_IBL_S02', 'experiment_id': 'CNC_Archive_E02734', 'scan': 7}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3e-06, 'EchoTime': 0.00169, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.00886408, 0.997264, -0.0733949, 0.021191, -0.0735686, -0.997065], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 3, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 93.75, 'PercentSampling': 100, 'PhaseEncodingSteps': 238, 'PhaseResolution': 1, 'PixelBandwidth': 650, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p3', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7;NC1', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 240, 'RefLinesPE': 32, 'RepetitionTime': 2.53, 'SAR': 0.0370875, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_16ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p3 RMS', 'SeriesNumber': 7, 'ShimSetting': [-925, -7137, 3764, 139, -99, -646, -157, -21], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 240.054, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot', 'dataset': '<unset>'}, 'provenance': {'md5sum': 'c9652695037ed45b6a0df330b0c813c6', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}}
* run_id : <undefined>
* session_id : 231101IBLS02
* subject_id : 231101IBLS02
* task_id : <undefined>


Execution Outputs
-----------------


* out_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/mriqc/sub-231101IBLS02/ses-231101IBLS02/anat/sub-231101IBLS02_ses-231101IBLS02_T1w.json


Runtime info
------------


* duration : 0.004918
* hostname : scc-xk3
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231101IBLS02..ses-231101IBLS02..anat..sub-231101IBLS02_ses-231101IBLS02_T1w.nii.gz/datasink


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

