function single_power_response_ft_correctVSwrong(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite,iterN)
% IsBL2preDelay: use pre-trial interval or pre response as baseline; better named as IsBL2preResp
rng shuffle
IsCorretTrials = 0;% must be 0, use all trials
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile = fullfile(Dir.results,[subname,'_respPower_ft_CorrWrong',num2str(iterN),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
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
    for cond_i = 3%1:3
        tmpID = find(M.ss_num == ss(cond_i));
        condTrials = M(tmpID,:);

        sEEG =  pop_select(EEG,'trial',tmpID);

        bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
        timelimits = [-0.88 1.5];% pre-stimuli baseline
        [bEEG,indices] = pop_epoch(bEEG,stim1,timelimits);
        rmv = setdiff(1:sEEG.trials,indices);
        sEEG = pop_select(sEEG,'notrial',rmv);
        condTrials(rmv,:) = [];

        sEEG =  pop_select(sEEG,'time',[-0.5 inf]);% shorten segment
        origTrials = 1:sEEG.trials;
        [sEEG,indx] = pop_epoch(sEEG,keyMarker,[-2 1]);% re-epoch to reponse
        rmvd = setdiff(origTrials,indx);% get removed epochs
        condTrials(rmvd,:) = [];

        bEEG = pop_select(bEEG,'notrial',rmvd);% remove from baselin data
        if sEEG.trials~=bEEG.trials% check trials N
            error('trials missing')
        end

        %%
        eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');
        beeg = eeglab2fieldtrip(bEEG,'preprocessing','none');

        if IsdePhase
            avgTrl = ft_timelockanalysis([], eeg);
            eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

            avgTrl = ft_timelockanalysis([], beeg);
            beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
        end

        wrongTrlN = sum(condTrials.button_resp_corr == 0);
        correctTrlN = sum(condTrials.button_resp_corr == 1);
        wrongID = find(condTrials.button_resp_corr ==0);
        correctID = find(condTrials.button_resp_corr ==1);
        %%
        c =1; % correct trials should be more than incorrect trials, therefore subsampling correct trials

        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.out = 'pow';
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.pad = 20;
        cfg.keeptrials  = 'yes';% dB normalization based on trial averaged baseline
        cfg.trials = correctID;
        tmp_eeg = ft_freqanalysis(cfg,eeg);

        if IsBL2preDelay==0 %baseline to pre trial
            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            cfg.keeptrials  = 'no';% dB normalization based on trial averaged baseline
            tmp_beeg = ft_freqanalysis(cfg,beeg);
        end

        if IsBL2preDelay
            cfg = [];
            cfg.baseline     = [-0.4 -0.1];% pre-response -0.4~-0.1s as baseline
            cfg.baselinetype = 'db';
            tmp_eeg = ft_freqbaseline(cfg,tmp_eeg);
        else
            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,tmp_beeg);

            timedim = find(size(bl.powspctrm)==length(bl.time));
            clear avgBL
            avgBL(1,:,:) = mean(bl.powspctrm,timedim);% add trial dimension

            tmp_eeg.powspctrm = 10*log10(bsxfun(@rdivide,tmp_eeg.powspctrm,avgBL)); %  db = 10*log10(signal/baseline);
        end

        clear pow

        for iter_i = 1:iterN
            tmpID = randperm(correctTrlN,wrongTrlN);

            cfg = [];
            cfg.trials = tmpID;
            tmp = ft_freqdescriptives(cfg,tmp_eeg);

            pow(iter_i,:,:,:) = tmp.powspctrm;
        end
        tfDat{cond_i,c} = ft_freqdescriptives([],tmp_eeg);% get ft structure;
        tfDat{cond_i,c}.powspctrm = squeeze(mean(pow));
        %%
        c = 2;%incorrect

        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.out = 'pow';
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.pad = 20;
        cfg.trials = wrongID;
        tmp_eeg = ft_freqanalysis(cfg,eeg);

        if IsBL2preDelay==0 %baseline to pre trial
            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            tmp_beeg = ft_freqanalysis(cfg,beeg);
        end

        cfg = [];
        cfg.latency = [-0.5 inf];
        tmp_eeg = ft_selectdata(cfg,tmp_eeg);

        if IsBL2preDelay
            cfg = [];
            cfg.baseline       = [-0.4 -0.1];% pre-response -0.4~-0.1s as baseline
            cfg.baselinetype = 'db';
            tmp_eeg = ft_freqbaseline(cfg,tmp_eeg);
        else
            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,tmp_beeg);

            timedim = find(size(bl.powspctrm)==length(bl.time));
            tmp_eeg.powspctrm = 10*log10(bsxfun(@rdivide,tmp_eeg.powspctrm,mean(bl.powspctrm,timedim))); %  db = 10*log10(signal/baseline);
        end

        tfDat{cond_i,c} = ft_freqdescriptives([],tmp_eeg);% get ft structure;

    end
    save(outputFile,'tfDat','-v7.3')
end
