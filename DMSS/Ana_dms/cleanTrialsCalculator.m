clear
load('subs.mat');

trialN = nan(height(subs),1);
origN = nan(height(subs),1);
parfor sn = 1:height(subs)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
    if isfile(set_name) && subs.exclude(sn)==0
        EEG = pop_loadset('filename',set_name);

        trialN(sn) = EEG.trials;
        set_name = fullfile(Dir.prepro,[subname,'_epo.set']);
        EEG = pop_loadset('filename',set_name);
        origN(sn) = EEG.trials;
    end
end

subs.cleanTrial = trialN;
subs.oriN = origN;
subs.cleanRatio = subs.cleanTrial./subs.oriN;
save(fullfile(Dir.ana,'subs.mat'),'subs','Dir')


datastats(subs.cleanRatio(subs.group==1 & subs.exclude==0 & subs.rawEEG==1))

datastats(subs.cleanRatio(subs.group==2 & subs.exclude==0 & subs.rawEEG==1))
