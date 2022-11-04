function single_variablity_delay_ft(sn,IsLap,IsdePhase,IsCorretTrials)
load('subs.mat');
if subs.rawEEG(sn)==1

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_delay.set']);
    EEG = pop_loadset('filename',set_name);% baseline corrected to pre-stim
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

        if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
            stim1 = {'S 11','S 21', 'S 31'};
        else
            stim1 = {'11','21','31'};
        end

        bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
        timelimits = [-0.85 3];% pre-stimuli baseline
        bEEG = pop_epoch(bEEG,stim1,timelimits);

        if sEEG.trials~=bEEG.trials
            error('trials missing')
        end

        eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');
        beeg = eeglab2fieldtrip(bEEG,'preprocessing','none');

        if IsdePhase
            avgTrl = ft_timelockanalysis([], eeg);
            eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

            avgTrl = ft_timelockanalysis([], beeg);
            beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
        end

        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.out = 'pow';
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.keeptrials  = 'yes';
        eeg = ft_freqanalysis(cfg,eeg);

        if IsBL2preDelay==0 %baseline to pre stim
            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            beeg = ft_freqanalysis(cfg,beeg);
        end

        cfg = [];
        cfg.latency = [-0.4 inf];
        eeg = ft_selectdata(cfg,eeg);% select time of interest

        if IsBL2preDelay
            cfg = [];
            cfg.latency = [-0.4 -0.1];% select baseline as preDelay
            bl = ft_selectdata(cfg,eeg);

            nv = squeeze(log(var(eeg.powspctrm,0,1)./mean(var(bl.powspctrm,0,1),4)));%chan*freq*time

        else %baseline to pre stim

            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,beeg);% select baseline as preStim

            nv = squeeze(log(var(eeg.powspctrm,0,1)./mean(var(bl.powspctrm,0,1),4)));%chan*freq*time
        end

        eeg = ft_freqdescriptives([],eeg);% get ft structure;
        eeg.powspctrm = nv;

        tfDat{cond_i} = eeg;
    end
    %%
    txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
    save(fullfile(Dir.results,[subname,'_var_delay_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'tfDat','-v7.3')

end
