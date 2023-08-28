clear
close all
load('subs.mat');

for sn = 1:height(subs)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_epo.set']);
    if isfile(set_name) && subs.exclude(sn)==0
        EEG = pop_loadset('filename',set_name);

        %%

        EEG = pop_select(EEG,'time',[-2.2 EEG.xmax]);% re-epoch into shorter segments

        csvFile = fullfile(Dir.beha,subs.csvFile{sn});
        M = readtable(csvFile);
        M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});

        M(isnan(M.ss_num),:) = [];
        if sum(isnan(M.button_resp_rt)) >0
            M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
        end

        tmp = find(contains({EEG.event.type},'T'));
        if length(tmp) == EEG.trials
            cleanTrials = {EEG.event(tmp).type};
            cleanTrials = cellfun(@(x)(str2num(x(2:end))),cleanTrials,'UniformOutput',false);
            cleanTrials = cell2mat(cleanTrials);
        end
        M = M(cleanTrials,:);

        checkMarkerBeha(EEG,M)

    end
end