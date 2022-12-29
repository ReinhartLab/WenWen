Intro:
	task following "Frontal Oscillatory Dynamics Predict Feedback Learning and Action Adjustment"
	data is offline located at D:\rob\timeest

	experiment code by John? in github ReinhartLab/Nona/Experiments/Reinhart Lab Experiments/Time Estimation/
	analysis code by Wen in github ReinhartLab/WenWen/TimeEst

	
baseline = -0.5~-0.1 precue, ROI = post feedback 0.1~0.5, FCz-F6, PLV

marker codes:

	Exp_beg = 1; %begin main experiment
	Exp_end = 2; %experiment end
	Test_Beg = 3;
	Test_End = 4; %begin entire test

	Block_New = 7; %new block = new time selected - blocks are essentially the "j" variable
	Block_Feedback = 8; %first time selected
	Block_NonFeedback = 9;

	Block_Learned = 33; %denotes whether the preceding block was learned (i.e. when success counter breaks first)
	Block_UnLearned = 34; %if preceding blocked was not learned (i.e. when error counter breaks first)

	Trial_start = 10; %start of each trial  (600 trials total). essentially the "i" variable
	Trial_end = 11; %end of each trial

	Fixation_Onset = 35; %Fixation onset - fixation cross
	Fixation_Offset = 36; %fixation offset - fixation cross

	Stim_Onset = 12; %stim onset - cue box
	Stim_Offset = 13; %stim offset - cue box disappears at same time as fixation cross

	Resp_all  = 17; % General Response Marker (when a participant responds)
	Resp_none = 18; % no response

	Resp_early = 19; %early response
	Resp_late = 20;    %late response
	Resp_correct = 21; %correct response

	Onset_fdb = 27; % General Feedback onset
	Fdb_correct   = 28; % Feedback Correct
	Fdb_early  = 29; % Feedback too early
	Fdb_late  = 30; % Feedback too late
	Offset_fdb = 31;%offset feedback
	fdb_incorrect = 37;
	Break_onset = 32; %30 second break onset

trialwise marker codes:
	%may have
		block start:5
	
	%must have 
			trial start: 10
			fixation: 35
			cue onset: 12
			action: 17-yes; 18-no 
			cue offset: 13
			respType:19-early; 20-later; 21-correct
			feedback onset:29-early; 30-late; 28-correct
			trial end:11

	%may have
		block: 7/8/9
		preceding block: 33-learned; 34-not learned
		block end:6

preprocessed by Wen and Doug:
	sn = 8,  S067, 350 trials, excluded from analysis
	sn = 10, S069, trial 100~400, 500uV, excluded from analysis
	sn = 14, S075, raw EEG.event(2531:2537) = [], eeglab 'boundary'

status logs:
	2022-10-20-wen: run loadRaw for everyone
	2022-10-20-doug: run intpChans for everyone
	2022-10-24-wen: fixed intpChans' error
	2022-10-27-wen: changed preprocessing order, epoch then interpolate chans
	2022-11-3-doug: doICA for sn=1
	2022-11-5-doug: postICA for sn=1
	2022-11-5-doug: doICA for sn=2:24
	2022-11-20-doug: postICA for sn=2:24
	2022-11-21-wen: check ICA for sn= 1 2 4:7 9 11 12 15:24
	2022-11-23-wen: check ICA for sn= 13 14, rerun preprocess sn=3
	2022-12-12-wen: sanity check of all data, sn13 noisy

notes for sn data (preICA):
	sn=3, first 60 trials are very noisy, the rest is noisy
	sn=6, very noisy
	sn=9, a lot of alpha waves; RHEOG & LHEOG broken
	sn=12, BVEOG broken
	sn=13, beautiful data, FCz seems to have opposite pattern.
	sn=15, removed multiple large distortions in the data; 434 trials left
	sn=16, very noisy
	sn=21, frontal channels are very noisy
	sn=23, a lot of alpha waves

notes for sn data (postICA):
	sn=3, pretty noisy
	sn=6, a bit noisy, removed many ICs
	sn=9, too many trials rejected (0.585 remaining), alpha
	sn=13, strange topology ICA components
	sn=15, removed many ICs, 0.66 trials remaining
	sn=16, a bit noisy (alpha waves)
	sn=20, a bit noisy
	sn=21, very noisy + alpha waves
	sn=22, possibly some eye artifacts left that are characterized as brain

