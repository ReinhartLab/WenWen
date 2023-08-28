clear
load('subs.mat');

trialN = nan(height(subs),1);
parfor sn = 1:height(subs)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
    if isfile(set_name) && subs.exclude(sn)==0
        EEG = pop_loadset('filename',set_name);
        trialN(sn) = EEG.trials;
    end
end

subs.cleanTrial = trialN
save(fullfile(Dir.ana,'subs.mat'),'subs','Dir')
