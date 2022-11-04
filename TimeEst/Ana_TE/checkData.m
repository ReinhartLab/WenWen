function checkData(sn)
threshV = 100;
locthresh = 8;
globthresh = 8;
zthresh = 8;

load('subs.mat');

subname = subs.name{sn};
set_name = fullfile(Dir.prepro,[subname,'_postICA.set']);
EEG = pop_loadset('filename',set_name);

pop_eegplot(EEG,1,1,1);

%% rejection
sEEG = pop_select(EEG,'time',[-0.5 4],'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

[~, rejT] = pop_eegthresh( sEEG, 1, 1:sEEG.nbchan, -threshV, threshV, sEEG.xmin, sEEG.xmax,0,0);
[~, ~, ~, rejK] = pop_rejkurt_ww( sEEG, 1, 1:sEEG.nbchan, locthresh,globthresh, 'reject',0,'plotflag',0);
[~, ~, ~, rejP] = pop_jointprob_ww( sEEG, 1, 1:sEEG.nbchan, locthresh, globthresh,1, 'reject',0,'plotflag',0);

tmp = abs(zscore(sEEG.data,0,3))>zthresh;% along trial dimension
rejZ = squeeze(sum(sum(tmp),2))>0;

tmp = abs(zscore(sEEG.data,0,2))>zthresh;% along time dimension
rejZ2 = squeeze(sum(sum(tmp),2))>0;

%% summary
trialN = 600;
rej_union = union(rejT,[find(rejK) find(rejP) find(rejZ)' find(rejZ2)']);
fprintf(['\n------trials summary-----\n' ...
    'Rejecting %d trials of %d\n' ...
    'Remaining %d of %d (%.3f)\n'], ...
    length(rej_union),EEG.trials,EEG.trials-length(rej_union),trialN,(EEG.trials-length(rej_union))/trialN);

%% remove

EEG = pop_rejepoch( EEG, rej_union,0);
pop_eegplot(EEG,1,1,1);

%% save

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
pop_saveset(EEG,'filename',set_name);

end
