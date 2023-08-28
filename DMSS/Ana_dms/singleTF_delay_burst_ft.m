function singleTF_delay_burst_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_delay_ft_burst',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

if IsOverwrite==0 && isfile(outputFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
if isfile(set_name)
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
    if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
        stim1 = {'S 11','S 21', 'S 31'};
    else
        stim1 = {'11','21','31'};
    end

    clear tfDat
    ss = [1 2 4];% set size
    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG =  pop_select(EEG,'trial',find(tmpID));
        eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');

        bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
        timelimits = [-0.88 3];% pre-stimuli baseline
        bEEG = pop_epoch(bEEG,stim1,timelimits);

        if sEEG.trials~=bEEG.trials
            error('trials missing')
        end

        beeg = eeglab2fieldtrip(bEEG,'preprocessing','none');

        if IsdePhase
            avgTrl = ft_timelockanalysis([], eeg);
            eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

            avgTrl = ft_timelockanalysis([], beeg);
            beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
        end

        cfg               = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi        = 1:40;
        cfg.pad           = 20;% in s
        cfg.tapsmofrq     = 1;
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.6;   % length of time window = 0.5 sec
        cfg.keeptrials  = 'yes';
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        eeg = ft_freqanalysis(cfg, eeg);

        if IsBL2preDelay
            %             cfg = [];
            %             cfg.baseline       = [-0.4 -0.1];% pre-delay -0.4~-0.1s as baseline
            %             cfg.baselinetype   = 'db';
            %             eeg = ft_freqbaseline(cfg, eeg);

        else % to pre trial

            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            beeg = ft_freqanalysis(cfg,beeg);

            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,beeg);

            tmp_data = squeeze(mean(bl.powspctrm,4)); % average time

            eeg.powspctrm = bsxfun(@(x,y)((x-mean(y,'omitnan'))/std(y,0,1,'omitnan')),eeg.powspctrm,tmp_data);% n SD: (x-bl)/sd

        end

        tfDat{cond_i} = eeg;
    end
    %%
    save(outputFile,'tfDat','-v7.3')

end
