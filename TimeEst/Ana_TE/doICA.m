function doICA(sn)

load('subs.mat')

subname = subs.name{sn};
set_name = fullfile(Dir.prepro,[subname,'_del.set']);

EEG = pop_loadset('filename',set_name);

%% https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html

% Start with an unfiltered (or minimally filtered) dataset (dataset 1)
% Filter the data at 1Hz or 2Hz to obtain dataset 2
% Run ICA on dataset 2
% Apply the resulting ICA weights to dataset 1
badChanFile = fullfile(Dir.ana,'BadChans',[subname,'_badchans.mat']);

if isfile(badChanFile)
    load(badChanFile);
    badchanN = length(badChan.interped);
else
    badchanN = 0;
end

if sn == 13 || badchanN>=6

    EEGica = pop_eegfiltnew(EEG,2);
    EEGica = pop_runica(EEGica,'icatype','runica','extended',1,'stop', 1E-7,'pca',50);%
    EEG.icaweights = EEGica.icaweights;
    EEG.icaact = EEGica.icaact;
    EEG.icachansind = EEGica.icachansind;
    EEG.icasphere = EEGica.icasphere;
    EEG.icawinv = EEGica.icawinv;

else

    EEG = pop_runica(EEG,'icatype','runica','extended',1,'stop', 1E-7);
end

%%
new_name = [subname,'_ICA.set'];
pop_saveset(EEG,'filename',new_name,'filepath',Dir.prepro);
end

