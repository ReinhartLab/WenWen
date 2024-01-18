Node: anatMRIQC (ComputeIQMs (datasink (bids)
=============================================


 Hierarchy : mriqc_wf.anatMRIQC.ComputeIQMs.datasink
 Exec ID : datasink.a0


Original Inputs
---------------


* _outputs : {'qi_2': 0.004792821411350609}
* acq_id : <undefined>
* dataset : <unset>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/sub-231117IBLS04/ses-231117IBLS04/anat/sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz
* metadata : {'AcquisitionMatrixPE': 230, 'AcquisitionNumber': 1, 'AcquisitionTime': '14:13:7.040000', 'BaseResolution': 230, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231117_IBL_S04', 'subject_id': 'CNC_Archive_S02567', 'experiment': '231117_IBL_S04', 'experiment_id': 'CNC_Archive_E02780', 'scan': 8}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3.3e-06, 'EchoTime': 0.00167, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.0675592, 0.996116, 0.0564728, -0.0363975, 0.0590251, -0.997593], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 2, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 100, 'PercentSampling': 100, 'PhaseEncodingSteps': 230, 'PhaseResolution': 1, 'PixelBandwidth': 660, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p2', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 230, 'RefLinesPE': 32, 'RepetitionTime': 2.2, 'SAR': 0.0455766, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_176ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p2 RMS', 'SeriesNumber': 8, 'ShimSetting': [-924, -7109, 3732, 304, -8, -644, -89, -10], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 216.42, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot'}
* modality : T1w
* out_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/mriqc
* provenance : {'md5sum': 'abc8cd9ff6a6eb2132e0d4ec8be815fc', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}
* rec_id : <undefined>
* root : {'summary_csf_mean': 276.725153247976, 'summary_csf_stdv': 111.6529666516892, 'summary_csf_median': 261.73464741557837, 'summary_csf_mad': 93.71215172783202, 'summary_csf_p95': 486.11225739121437, 'summary_csf_p05': 106.4232488758862, 'summary_csf_k': 2.0493367802475806, 'summary_csf_n': 38142.0, 'summary_gm_mean': 732.8021734237298, 'summary_gm_stdv': 56.47424245101941, 'summary_gm_median': 733.114552013576, 'summary_gm_mad': 54.29762131875481, 'summary_gm_p95': 825.9560670033095, 'summary_gm_p05': 639.2636493183672, 'summary_gm_k': 0.23096545556355919, 'summary_gm_n': 37570.0, 'summary_wm_mean': 1000.9927181829861, 'summary_wm_stdv': 39.03453666205313, 'summary_wm_median': 1000.027676820755, 'summary_wm_mad': 36.58161103592479, 'summary_wm_p95': 1067.541729658842, 'summary_wm_p05': 938.7721046805382, 'summary_wm_k': 0.47688328698068894, 'summary_wm_n': 172811.0, 'summary_bg_mean': 8.396000727183301, 'summary_bg_stdv': 19.51401914605082, 'summary_bg_median': 4.707706719636917, 'summary_bg_mad': 6.9796564266074235, 'summary_bg_p95': 26.75085112452507, 'summary_bg_p05': 0.0, 'summary_bg_k': 1093.7797767320044, 'summary_bg_n': 2606239.0, 'snr_csf': 2.344149234779899, 'snr_wm': 25.618974091072396, 'snr_gm': 12.981224067524177, 'snr_total': 13.648115797792158, 'snrd_csf': 24.567383289732792, 'snrd_wm': 93.86630115419041, 'snrd_gm': 68.81284679899922, 'snrd_total': 62.41551041430747, 'cnr': 3.7397934670471447, 'fber': 5833.6967668960315, 'efc': 0.5310809478016508, 'wm2max': 0.68122995661196, 'qi_1': 0.00015013877163003387, 'cjv': 0.3404824413199305, 'fwhm_x': 4.21469, 'fwhm_y': 4.06536, 'fwhm_z': 3.90834, 'fwhm_avg': 4.062796666666666, 'icvs_csf': 0.2051232616003943, 'icvs_gm': 0.4538435819557907, 'icvs_wm': 0.341033156443815, 'rpve_csf': 22.77205889597147, 'rpve_gm': 10.982222902787033, 'rpve_wm': 16.202051132838108, 'size_x': 176, 'size_y': 230, 'size_z': 230, 'spacing_x': 1.0, 'spacing_y': 1.0, 'spacing_z': 1.0, 'inu_range': 0.1991374075412749, 'inu_med': 0.4928629696369171, 'tpm_overlap_csf': 0.14293616198343354, 'tpm_overlap_gm': 0.47367844201119064, 'tpm_overlap_wm': 0.5027335095491199, 'qi_2': 0.004792821411350609, 'bids_meta': {'subject_id': '231117IBLS04', 'session_id': '231117IBLS04', 'modality': 'T1w', 'AcquisitionMatrixPE': 230, 'AcquisitionNumber': 1, 'AcquisitionTime': '14:13:7.040000', 'BaseResolution': 230, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231117_IBL_S04', 'subject_id': 'CNC_Archive_S02567', 'experiment': '231117_IBL_S04', 'experiment_id': 'CNC_Archive_E02780', 'scan': 8}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3.3e-06, 'EchoTime': 0.00167, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.0675592, 0.996116, 0.0564728, -0.0363975, 0.0590251, -0.997593], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 2, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 100, 'PercentSampling': 100, 'PhaseEncodingSteps': 230, 'PhaseResolution': 1, 'PixelBandwidth': 660, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p2', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 230, 'RefLinesPE': 32, 'RepetitionTime': 2.2, 'SAR': 0.0455766, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_176ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p2 RMS', 'SeriesNumber': 8, 'ShimSetting': [-924, -7109, 3732, 304, -8, -644, -89, -10], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 216.42, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot', 'dataset': '<unset>'}, 'provenance': {'md5sum': 'abc8cd9ff6a6eb2132e0d4ec8be815fc', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}}
* run_id : <undefined>
* session_id : 231117IBLS04
* subject_id : 231117IBLS04
* task_id : <undefined>


Execution Inputs
----------------


* _outputs : {'qi_2': 0.004792821411350609}
* acq_id : <undefined>
* dataset : <unset>
* in_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/sub-231117IBLS04/ses-231117IBLS04/anat/sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz
* metadata : {'AcquisitionMatrixPE': 230, 'AcquisitionNumber': 1, 'AcquisitionTime': '14:13:7.040000', 'BaseResolution': 230, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231117_IBL_S04', 'subject_id': 'CNC_Archive_S02567', 'experiment': '231117_IBL_S04', 'experiment_id': 'CNC_Archive_E02780', 'scan': 8}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3.3e-06, 'EchoTime': 0.00167, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.0675592, 0.996116, 0.0564728, -0.0363975, 0.0590251, -0.997593], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 2, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 100, 'PercentSampling': 100, 'PhaseEncodingSteps': 230, 'PhaseResolution': 1, 'PixelBandwidth': 660, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p2', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 230, 'RefLinesPE': 32, 'RepetitionTime': 2.2, 'SAR': 0.0455766, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_176ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p2 RMS', 'SeriesNumber': 8, 'ShimSetting': [-924, -7109, 3732, 304, -8, -644, -89, -10], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 216.42, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot'}
* modality : T1w
* out_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/mriqc
* provenance : {'md5sum': 'abc8cd9ff6a6eb2132e0d4ec8be815fc', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}
* rec_id : <undefined>
* root : {'summary_csf_mean': 276.725153247976, 'summary_csf_stdv': 111.6529666516892, 'summary_csf_median': 261.73464741557837, 'summary_csf_mad': 93.71215172783202, 'summary_csf_p95': 486.11225739121437, 'summary_csf_p05': 106.4232488758862, 'summary_csf_k': 2.0493367802475806, 'summary_csf_n': 38142.0, 'summary_gm_mean': 732.8021734237298, 'summary_gm_stdv': 56.47424245101941, 'summary_gm_median': 733.114552013576, 'summary_gm_mad': 54.29762131875481, 'summary_gm_p95': 825.9560670033095, 'summary_gm_p05': 639.2636493183672, 'summary_gm_k': 0.23096545556355919, 'summary_gm_n': 37570.0, 'summary_wm_mean': 1000.9927181829861, 'summary_wm_stdv': 39.03453666205313, 'summary_wm_median': 1000.027676820755, 'summary_wm_mad': 36.58161103592479, 'summary_wm_p95': 1067.541729658842, 'summary_wm_p05': 938.7721046805382, 'summary_wm_k': 0.47688328698068894, 'summary_wm_n': 172811.0, 'summary_bg_mean': 8.396000727183301, 'summary_bg_stdv': 19.51401914605082, 'summary_bg_median': 4.707706719636917, 'summary_bg_mad': 6.9796564266074235, 'summary_bg_p95': 26.75085112452507, 'summary_bg_p05': 0.0, 'summary_bg_k': 1093.7797767320044, 'summary_bg_n': 2606239.0, 'snr_csf': 2.344149234779899, 'snr_wm': 25.618974091072396, 'snr_gm': 12.981224067524177, 'snr_total': 13.648115797792158, 'snrd_csf': 24.567383289732792, 'snrd_wm': 93.86630115419041, 'snrd_gm': 68.81284679899922, 'snrd_total': 62.41551041430747, 'cnr': 3.7397934670471447, 'fber': 5833.6967668960315, 'efc': 0.5310809478016508, 'wm2max': 0.68122995661196, 'qi_1': 0.00015013877163003387, 'cjv': 0.3404824413199305, 'fwhm_x': 4.21469, 'fwhm_y': 4.06536, 'fwhm_z': 3.90834, 'fwhm_avg': 4.062796666666666, 'icvs_csf': 0.2051232616003943, 'icvs_gm': 0.4538435819557907, 'icvs_wm': 0.341033156443815, 'rpve_csf': 22.77205889597147, 'rpve_gm': 10.982222902787033, 'rpve_wm': 16.202051132838108, 'size_x': 176, 'size_y': 230, 'size_z': 230, 'spacing_x': 1.0, 'spacing_y': 1.0, 'spacing_z': 1.0, 'inu_range': 0.1991374075412749, 'inu_med': 0.4928629696369171, 'tpm_overlap_csf': 0.14293616198343354, 'tpm_overlap_gm': 0.47367844201119064, 'tpm_overlap_wm': 0.5027335095491199, 'qi_2': 0.004792821411350609, 'bids_meta': {'subject_id': '231117IBLS04', 'session_id': '231117IBLS04', 'modality': 'T1w', 'AcquisitionMatrixPE': 230, 'AcquisitionNumber': 1, 'AcquisitionTime': '14:13:7.040000', 'BaseResolution': 230, 'BodyPartExamined': 'BRAIN', 'ConsistencyInfo': 'N4_VE11C_LATEST_20160120', 'ConversionSoftware': 'dcm2niix', 'ConversionSoftwareVersion': 'v1.0.20210317', 'DataSource': {'application/x-xnat': {'url': 'https://xnat.bu.edu', 'project': 'Reinhart_IBL', 'subject': '231117_IBL_S04', 'subject_id': 'CNC_Archive_S02567', 'experiment': '231117_IBL_S04', 'experiment_id': 'CNC_Archive_E02780', 'scan': 8}}, 'DeviceSerialNumber': '166024', 'DwellTime': 3.3e-06, 'EchoTime': 0.00167, 'FlipAngle': 7, 'ImageOrientationPatientDICOM': [0.0675592, 0.996116, 0.0564728, -0.0363975, 0.0590251, -0.997593], 'ImageType': ['ORIGINAL', 'PRIMARY', 'OTHER', 'ND', 'NORM', 'MEAN'], 'ImagingFrequency': 123.249, 'InPlanePhaseEncodingDirectionDICOM': 'ROW', 'InstitutionAddress': 'Commonwealth Ave. 610,Boston,Massachusetts,US,02215', 'InstitutionName': 'BU', 'InstitutionalDepartmentName': 'Research', 'InversionTime': 1.1, 'MRAcquisitionType': '3D', 'MagneticFieldStrength': 3, 'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Prisma', 'Modality': 'MR', 'NumberOfAverages': 4, 'ParallelReductionFactorInPlane': 2, 'PartialFourier': 1, 'PatientPosition': 'HFS', 'PercentPhaseFOV': 100, 'PercentSampling': 100, 'PhaseEncodingSteps': 230, 'PhaseResolution': 1, 'PixelBandwidth': 660, 'ProcedureStepDescription': 'INVESTIGATORS^Reinhart', 'ProtocolName': 'T1_MEMPRAGE_1.0mm_p2', 'PulseSequenceDetails': '%CustomerSeq%\\tfl_mgh_epinav_ABCD', 'ReceiveCoilActiveElements': 'HC1-7', 'ReceiveCoilName': 'HeadNeck_64', 'ReconMatrixPE': 230, 'RefLinesPE': 32, 'RepetitionTime': 2.2, 'SAR': 0.0455766, 'ScanOptions': 'IR', 'ScanningSequence': 'GR\\IR', 'SequenceName': 'tfl3d4_176ns', 'SequenceVariant': 'SK\\SP\\MP', 'SeriesDescription': 'T1_MEMPRAGE_1.0mm_p2 RMS', 'SeriesNumber': 8, 'ShimSetting': [-924, -7109, 3732, 304, -8, -644, -89, -10], 'SliceThickness': 1, 'SoftwareVersions': 'syngo MR E11', 'StationName': 'AWP166024-4CA74C', 'TxRefAmp': 216.42, 'WipMemBlock': 'Prisma_epi_moco_navigator_ABCD_tfl.prot', 'dataset': '<unset>'}, 'provenance': {'md5sum': 'abc8cd9ff6a6eb2132e0d4ec8be815fc', 'version': '22.0.6', 'software': 'mriqc', 'webapi_url': 'https://mriqc.nimh.nih.gov/api/v1', 'webapi_port': None, 'settings': {'testing': False}, 'warnings': {'small_air_mask': False, 'large_rot_frame': False}}}
* run_id : <undefined>
* session_id : 231117IBLS04
* subject_id : 231117IBLS04
* task_id : <undefined>


Execution Outputs
-----------------


* out_file : /projectnb/viscog01/Wen/IBL_TI_fMRI/BIDS/derivatives/mriqc/sub-231117IBLS04/ses-231117IBLS04/anat/sub-231117IBLS04_ses-231117IBLS04_T1w.json


Runtime info
------------


* duration : 0.00428
* hostname : scc-xh4
* prev_wd : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code
* working_dir : /projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/work/mriqc_wf/anatMRIQC/ComputeIQMs/_in_file_..projectnb..viscog01..Wen..IBL_TI_fMRI..BIDS..sub-231117IBLS04..ses-231117IBLS04..anat..sub-231117IBLS04_ses-231117IBLS04_T1w.nii.gz/datasink


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
* TMPDIR : /scratch/2657782.1.onrcc-m256
* USER : wenwen

