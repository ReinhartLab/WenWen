function single_power_DelayResponse(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
%use pre-trial baseline for maintenance
% pre-response baseline for deletion
IsBL2preDelay = 1;% enforced

load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_delayresp_trialwise',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
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
    M.trl = [1:height(M)]';
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

    if IsCorretTrials
        tmpIdx = M.button_resp_corr ==1;
        M = M(tmpIdx,:);
        EEG = pop_select(EEG,'trial',find(tmpIdx));
    end

    sEEG = EEG;
    bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
    timelimits = [-0.88 1.5];% pre-stimuli baseline
    [bEEG,indices] = pop_epoch(bEEG,stim1,timelimits);
    rmv = setdiff(1:sEEG.trials,indices);
    sEEG = pop_select(sEEG,'notrial',rmv);
    M(rmv,:) = [];

    sEEG =  pop_select(sEEG,'time',[-0.5 inf]);% shorten segment
    origTrials = 1:sEEG.trials;
    [rEEG,indx] = pop_epoch(sEEG,keyMarker,[-3 0.8]);% re-epoch to reponse
    rmvd = setdiff(origTrials,indx);% get removed epochs
    M(rmvd,:) = [];

    bEEG = pop_select(bEEG,'notrial',rmvd);% remove from baseline data
    sEEG  = pop_select(sEEG,'notrial',rmvd);% remove from baseline data

    if sEEG.trials~=bEEG.trials% check trials N
        error('trials missing')
    end

    eeg = eeglab2fieldtrip(sEEG,'preprocessing','none'); % delay
    reeg = eeglab2fieldtrip(rEEG,'preprocessing','none'); % response
    beeg = eeglab2fieldtrip(bEEG,'preprocessing','none');%pre trial baseline

    if IsdePhase
        avgTrl = ft_timelockanalysis([], eeg);
        eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

        avgTrl = ft_timelockanalysis([], reeg);
        reeg.trial = cellfun(@(x)x-avgTrl.avg,reeg.trial,'UniformOutput',false);

        avgTrl = ft_timelockanalysis([], beeg);
        beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
    end
%%
    cfg = [];
    cfg.method = 'wavelet';
    cfg.foi = 1:40;
    cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
    cfg.out = 'pow';
    cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
    cfg.keeptrials  = 'yes';
    cfg.pad = 20;
    eeg = ft_freqanalysis(cfg,eeg);

    cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
    beeg = ft_freqanalysis(cfg,beeg);

    cfg.toi = rEEG.xmin:0.05:rEEG.xmax; % every 50ms
    reeg = ft_freqanalysis(cfg,reeg);

    cfg = [];
    cfg.latency = [-0.5 3.2];
    tfDat.PowN = ft_selectdata(cfg,eeg); % maintenance power

    %baseline to pre trial
    cfg = [];
    cfg.latency = [-0.4 -0.1];
    bl = ft_selectdata(cfg,beeg);
    % db-baseline trialwise power data
    clear avgBL
    timedimB = find(size(bl.powspctrm)==length(bl.time));
    avgBL = mean(mean(bl.powspctrm,1),timedimB,"omitnan");%  trial (1)*chan*freq
    tfDat.PowN.powspctrm = 10*log10(bsxfun(@rdivide,tfDat.PowN.powspctrm,avgBL)); %  db = 10*log10(signal/baseline);

    % baseline to pre-response
    cfg = [];
    cfg.baseline     = [-0.4 -0.1];% pre-response -0.4~-0.1s as baseline
    cfg.baselinetype = 'db';
    tfDat.PowN_1 = ft_freqbaseline(cfg,reeg); % resp power

    cfg = [];
    cfg.latency = [-0.4 0.8];
    tfDat.PowN_1 = ft_selectdata(cfg,tfDat.PowN_1);

    save(outputFile,'M','tfDat','-v7.3')
end
