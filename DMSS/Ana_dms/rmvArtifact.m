clear
load('subs.mat');
close all

sn = 17;%change this 
%%
if subs.rawEEG(sn) 
    %% load the data
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    subname = subs.name{sn};

    set_name = fullfile(Dir.prepro,[subname,'_epo.set']);
    EEG = pop_loadset('filename',set_name);

    %% reject trials by visual inspection, repeat until remove all bad trials

    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % Store dataset
    pop_eegplot(EEG,1,1,1);% plot the EEG time series.
    pause
    %% modified pop_rejkurt to [EEG, locthresh, globthresh, rej,nrej,com]

    sEEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

    threshV = 600;
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

    EEG = pop_rejepoch(EEG, rmvID,0);

    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); % Store dataset
    pop_eegplot(EEG,1,1,1);
    pause
    %%  save

    set_name = fullfile(Dir.prepro,[subname,'_del.set']);
    pop_saveset(EEG,'filename',set_name);

end
