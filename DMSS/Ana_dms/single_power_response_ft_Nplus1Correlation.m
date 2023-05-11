function single_power_response_ft_Nplus1Correlation(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
% IsBL2preDelay: use pre-trial interval or pre response as baseline; better named as IsBL2preResp
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_resp_PowerTrialwise',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
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
    if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
        stim1 = {'S 11','S 21', 'S 31'};
        keyMarker = {'S 60'};
    else
        stim1 = {'11','21','31'};
        keyMarker = {'60'};
    end

    clear tfDat
    ss = [1 2 4];% set size
    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        condTrials = find(tmpID);
        sEEG =  pop_select(EEG,'trial',condTrials);

        bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
        timelimits = [-0.88 1.5];% pre-stimuli baseline
        [bEEG,indices] = pop_epoch(bEEG,stim1,timelimits);
        rmv = setdiff(1:sEEG.trials,indices);
        sEEG = pop_select(sEEG,'notrial',rmv);
        condTrials(rmv) = [];

        sEEG =  pop_select(sEEG,'time',[-0.5 inf]);% shorten segment
        origTrials = 1:sEEG.trials;
        [sEEG,indx] = pop_epoch(sEEG,keyMarker,[-2 1]);% re-epoch to reponse
        rmvd = setdiff(origTrials,indx);% get removed epochs
        condTrials(rmvd) = [];

        bEEG = pop_select(bEEG,'notrial',rmvd);% remove from baselin data
        if sEEG.trials~=bEEG.trials% check trials N
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
        cfg.pad = 20;
        eeg = ft_freqanalysis(cfg,eeg);

        if IsBL2preDelay==0 %baseline to pre trial
            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            cfg.keeptrials  = 'no';% dB normalization based on trial averaged baseline
            beeg = ft_freqanalysis(cfg,beeg);
        end

        cfg = [];
        cfg.latency = [-1 inf];
        eeg = ft_selectdata(cfg,eeg);

        if IsBL2preDelay
            cfg = [];
            cfg.baseline       = [-0.4 -0.1];% pre-response -0.4~-0.1s as baseline
            cfg.baselinetype   = 'db';
            eeg = ft_freqbaseline(cfg, eeg);

        else %baseline to pre trial

            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,beeg);

            timedim = find(size(bl.powspctrm)==length(bl.time));
            clear avgBL
            avgBL(1,:,:) = mean(bl.powspctrm,timedim);% add trial dimension

            eeg.powspctrm = 10*log10(bsxfun(@rdivide,eeg.powspctrm,avgBL)); %  db = 10*log10(signal/baseline);
        end

        if condTrials(end)+1 > height(M) % no N+1 trial for the last trial
            cfg = [];
            cfg.trials = 1:(length(condTrials)-1);
            eeg = ft_selectdata(cfg,eeg);
            condTrials(end) = [];
        end

        tfDat{cond_i} = eeg;
        PostTrlBeha{cond_i} = M(condTrials+1,:);
    end
    save(outputFile,'tfDat','PostTrlBeha','-v7.3')
end
