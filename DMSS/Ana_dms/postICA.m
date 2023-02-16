clear;

load('subs.mat');
for sn = 44 %[20:21]%:height(subs)
    close all
    if subs.rawEEG(sn)
        subname = subs.name{sn};

        set_name = fullfile(Dir.prepro,[subname,'_ICA.set']);
        EEG = pop_loadset('filename',set_name);

        EEG = pop_iclabel(EEG,'default');
        pop_selectcomps(EEG, [1:min([40 size(EEG.icaweights,1)])]);
%         pop_selectcomps(EEG, [1:28]);
        pop_eegplot( EEG, 0, 1, 1);

        rejICA = fullfile(Dir.ana,'rmvICA',[subname,'_rmvID.mat']);
        if exist(rejICA)
            load(rejICA)
            disp(rmvd)
        end

        if isfile(fullfile(Dir.prepro,[subname,'_clean.set']))
            %   checkData(sn)
        end

        pause
        rmvd  = input('enter to-be-rmvd cmp ID with brackets[]:');

        EEG = pop_subcomp(EEG,rmvd,0); % change to 1 to confirm rejection

        set_name = fullfile(Dir.prepro,[subname,'_postICA.set']);
        pop_saveset(EEG,'filename',set_name);

        save(rejICA,'rmvd')
    end
end
