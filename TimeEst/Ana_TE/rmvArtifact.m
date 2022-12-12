clear
load('subs.mat');
close all

sn = 3;%change this

%% load the data

subname = subs.name{sn};

set_name = fullfile(Dir.prepro,[subname,'_pre.set']);
EEG = pop_loadset('filename',set_name);

pop_eegplot(EEG,1,1,1);% plot the EEG time series.

matname = fullfile(Dir.ana,'VisRemovalBeforeICA',[subname '_manualRmv.mat']);
if isfile(matname)
    load(matname)
    disp(rejID)
end

pause
rejID = input('enter to-be-rmvd trial ID with brackets[]:');
save(matname,'rejID')

%% modified pop_rejkurt to [EEG, locthresh, globthresh, rej,nrej,com]

sEEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels
sEEG = pop_rejepoch(sEEG, rejID,0);% remove inspected trials

threshV = 500;
[~, rejT] = pop_eegthresh( sEEG, 1, 1:sEEG.nbchan, -threshV, threshV, sEEG.xmin, sEEG.xmax,0,0);

locthresh = 10;
globthresh = 10;

[~, ~, ~, rejK] = pop_rejkurt_ww( sEEG, 1, 1:sEEG.nbchan, locthresh,globthresh, 'reject',0,'plotflag',0);
[~, ~, ~, rejP] = pop_jointprob_ww( sEEG, 1, 1:sEEG.nbchan, locthresh, globthresh,1, 'reject',0,'plotflag',0);

zthresh = 15;
tmp = abs(zscore(sEEG.data,0,3))>zthresh;% along trial dimension
rejZ = squeeze(sum(sum(tmp),2))>0;

tmp = abs(zscore(sEEG.data,0,2))>zthresh;% along time dimension
rejZ2 = squeeze(sum(sum(tmp),2))>0;

rmvID= union(rejT,[find(rejK) find(rejP)  find(rejZ)' find(rejZ2)']);

%% remove trials
EEG = pop_rejepoch(EEG, rejID,0);% reject visually identified epochs first
EEG = pop_rejepoch(EEG, rmvID,0);% then the automatic detected bad trials

pop_eegplot(EEG,1,1,1);

set_name = fullfile(Dir.prepro,[subname,'_del.set']);
pop_saveset(EEG,'filename',set_name);

