function singleTF_probe_plv_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
% IsBL2preDelay: use pre-trial interval or pre probe as baseline; better
% named as IsBL2preprobe
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_probe_ft_plv',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

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
        keyMarker = {'S 57'};%probe onset
    else
        stim1 = {'11','21','31'};
        keyMarker = {'57'};
    end

    clear tfDat
    ss = [1 2 4];% set size
    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG =  pop_select(EEG,'trial',find(tmpID));

        bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
        timelimits = [-0.88 1.5];% pre-stimuli baseline
        [bEEG,indices]= pop_epoch(bEEG,stim1,timelimits);
        rmv = setdiff(1:sEEG.trials,indices);
        sEEG = pop_select(sEEG,'notrial',rmv);

        sEEG =  pop_select(sEEG,'time',[-0.5 inf]);% shorten segment
        origTrials = 1:sEEG.trials;
        [sEEG,indx] = pop_epoch(sEEG,keyMarker,[-3 2]);% re-epoch to reponse
        rmvd = setdiff(origTrials,indx);% get removed epochs

        bEEG = pop_select(bEEG,'notrial',rmvd);% remove from baselin data
        if sEEG.trials~=bEEG.trials% check trials N
            error('trials missing')
        end

        beeg = eeglab2fieldtrip(bEEG,'preprocessing','none');
        eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');

        if IsdePhase
            avgTrl = ft_timelockanalysis([], eeg);
            eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

            avgTrl = ft_timelockanalysis([], beeg);
            beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
        end
        %% frequency decomposition
        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.output = 'fourier';
        cfg.pad = 20;
        freq = ft_freqanalysis(cfg,eeg);

        if IsBL2preDelay
            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bfreq = ft_selectdata(cfg,freq);

        else %baseline to pre trial

            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            bfreq = ft_freqanalysis(cfg,beeg);
        end

        %% plv
        cfg            = [];
        cfg.method     = 'plv';

        plv = ft_connectivityanalysis(cfg, freq);
        plvB = ft_connectivityanalysis(cfg, bfreq);

        cfg = [];
        cfg.avgovertime = 'yes';
        plvB = ft_selectdata(cfg,plvB);
        tfDat.plvB{cond_i} = plvB;% average time to save space

        cfg = [];
        cfg.latency = [-1 inf];% relative to probe
        tfDat.plv{cond_i} = ft_selectdata(cfg,plv);

%         %% wpli_debiased
%         cfg            = [];
%         cfg.method     = 'wpli';
% 
%         wpli = ft_connectivityanalysis(cfg, freq);
%         tfDat.wpliB{cond_i} = ft_connectivityanalysis(cfg, bfreq);
% 
%         cfg = [];
%         cfg.latency = [-1 inf];% tive to probe
%         tfDat.wpli{cond_i} =  ft_selectdata(cfg,wpli);
    end
    %%
    save(outputFile,'tfDat','-v7.3')
end
