%--------routine to set up environment
alias in terminal 
activate python vitual environment: va
intall packages (or skip): pip install xxx
save for TI: v_freeze
set for TI: vr_ti
end: deactivate


%-------STEP1: get data from XNAT to SCC in raw and BIDS format
ArcGetRaw_ww.sh 
ArcGetBIDS_ww.sh

%-------STEP2: run fmriprep on bids and perform QC
fmriprep_single_IBL_ww.qsub % need to change /projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/subj_single.txt, takes 12 hours+ depending on SCC status

mriqc_IBL.qsub

post_fmriprep_smooth_scaling.qsub # using fsl to run spatial smoothing and scaling; no need (wen 20240112)

fmriprep_multi_IBL_ww.qsub % run on all subjects in /projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/subj_all.txt, need to check requested time and cpus based on subjects number

%-------STEP3: behavioral analysis and generate timing files for fmri GLM
ppts_IBL.m % set up path and save subsinfo.mat
AnaBehaAll.m % analyzing behavior
AnaSingleSubj.m % analyzing individual data
genOnset.m % generate timing txt file for fsl 

%-------STEP4:localizer glm for FFA, OBA, MT
localizer_IBL.qsub

#-------------log 
20231103 wen 
	created ArcGetRaw_ww.sh 
20231104 wen 
	created ArcGetBIDS_ww.sh genOnset.m ppts_IBL.m
20231105 wen 
	created AnaBehaAll.m AnaSingleSubj.m
20231106 wen
	created fmriprep_single_IBL_ww.qsub mriqc_IBL.qsub
20231109 wen
post_fmriprep_smooth_scaling.qsub
