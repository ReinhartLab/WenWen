function single_SNR_delay(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

outputFile = fullfile(Dir.results,[subname,'_SNR_delay_',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

if IsOverwrite==0 && isfile(outputFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
if isfile(set_name)
    EEG = pop_loadset('filename',set_name);% baseline corrected to pre-trial
    EEG = pop_select(EEG,'nochannel',{'LHEOG','RHEOG','BVEOG','TVEOG'});
    if IsLap
        X = [EEG.chanlocs.X];
        Y = [EEG.chanlocs.Y];
        Z = [EEG.chanlocs.Z];
        EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
    end
    %% load behavior
    csvFile = fullfile(Dir.beha,subs.csvFile{sn});
    M = readtable(csvFile);
    M = M(:,["block_num","ss_num","button_resp_rt","button_resp_corr"]);
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
    if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
        stim1 = {'S 11','S 21', 'S 31'};
    else
        stim1 = {'11','21','31'};
    end
    SStime = {[-1.6 4],[-2.8 4],[-5.2 4]};% epoch length for different SS
    ss = [1 2 4];% set size
    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG = pop_select(EEG,'trial',find(tmpID),'time',SStime{cond_i}); % locked to delay onset

        if IsBL2preDelay==0 %baseline to pre trial
            [bEEG, indice]= pop_epoch(sEEG,stim1(cond_i),SStime{cond_i}(1)+[0 0.4]);
            sEEG = pop_select(sEEG,'trial',indice);
            sEEG.data = sEEG.data - repmat(mean(bEEG.data,2),[1 sEEG.pnts 1]);
        else %baseline to pre delay
            sEEG = pop_rmbase(sEEG,[-400 0]);
        end
        sEEG = pop_select(sEEG,'time',SStime{1}); %fixed data length for each set size

        SNR.broad(cond_i,:,:) = squeeze(mean(sEEG.data,3)./std(sEEG.data,0,3));% chan*time

        eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');

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
        cfg.pad = 20;

        eeg = ft_freqanalysis(cfg,eeg);

        SNR.freqdata(cond_i,:,:,:) = squeeze(mean(eeg.powspctrm)./std(eeg.powspctrm));% chan*freq*time
    end
    SNR.chans = sEEG.chanlocs;
    SNR.times = sEEG.times;
    SNR.ftoi = eeg.time;
    SNR.freq = eeg.freq;
    %%
    save(outputFile,'SNR','-v7.3')
end
