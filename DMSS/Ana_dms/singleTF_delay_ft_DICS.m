function    singleTF_delay_ft_DICS(sn,IsdePhase,IsCorretTrials,IsOverwrite)
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials'};
outputFile =fullfile(Dir.results,[subname,'_delay_ft',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']);

if IsOverwrite==0 && isfile(outputFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
if isfile(set_name)
    EEG = pop_loadset('filename',set_name);% baseline corrected to pre-stim
    EEG = pop_select(EEG,'nochannel',{'LHEOG','RHEOG','BVEOG','TVEOG'});
   
    %% load behavior
    csvFile = fullfile(Dir.beha,subs.csvFile{sn});
    M = readtable(csvFile);
    M = M(:,["block_num","ss_num","type","button_resp_rt","button_resp_corr"]);
    M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});

    M(isnan(M.ss_num),:) = [];
    if sum(isnan(M.button_resp_rt)) >0
        M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
    end

    tmp = find(contains({EEG.event.type},'T'));
    if length(tmp) ==EEG.trials
        goodTrials = {EEG.event(tmp).type};
        goodTrials = cellfun(@(x)(str2num(x(2:end))),goodTrials,'UniformOutput',false);
        goodTrials = cell2mat(goodTrials);
    end
    M = M(goodTrials,:);

    %%
    clear tfDat
    ss = [1 2 4];% set size
    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG =  pop_select(EEG,'trial',find(tmpID));
        eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');

        if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
            stim1 = {'S 11','S 21', 'S 31'};
        else
            stim1 = {'11','21','31'};
        end

      

        if IsdePhase
            avgTrl = ft_timelockanalysis([], eeg);
            eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);
        
        end

        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.out = 'pow';
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.keeptrials  = 'yes';
        eeg = ft_freqanalysis(cfg,eeg);

       
    end
    %%
    save(outputFile,'tfDat','-v7.3')

end
