function singleTF_response_ft_wavelet_EqualTrialSS(sn,iterN,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
% IsBL2preDelay: use pre-trial interval or pre response as baseline; better named as IsBL2preResp
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_resp_ft_equalizedTrials',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
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

    %% subsampling trials to match across set sizes

    [N,EDGES] = histcounts(M.ss_num,'BinWidth',1);
    %   N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1). The
    %   last bin will also include the right edge such that N(end) will count
    %   X(i) if EDGES(end-1) <= X(i) <= EDGES(end).

    ss = [1 2 4];% set size

    for iter = 1:iterN
        clear cond_idx tmp_rand rand_idx
        for i = 1:3
            cond_idx{i} = find(M.ss_num ==ss(i));
            tmp_rand{i} = randperm(N(i),min(N));
            rand_idx{i} =  cond_idx{i}(tmp_rand{i});
        end

        %%
        if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
            stim1 = {'S 11','S 21','S 31'};
            keyMarker = {'S 60'};
        else
            stim1 = {'11','21','31'};
            keyMarker = {'60'};
        end
        %%
        for cond_i = 1:3

            tmpID = rand_idx{cond_i};%get the randomized indices
            if IsCorretTrials ==1
                tmpID = intersect(find(M.button_resp_corr ==1),rand_idx{cond_i});
            end
            sEEG =  pop_select(EEG,'trial',tmpID);

            bEEG = pop_select(sEEG,'time',[-inf 2.5]);% re-epoch into shorter segments
            timelimits = [-0.88 1.5];% pre-stimuli baseline
            [bEEG,indices] = pop_epoch(bEEG,stim1,timelimits);
            rmv = setdiff(1:sEEG.trials,indices);
            sEEG = pop_select(sEEG,'notrial',rmv);

            sEEG =  pop_select(sEEG,'time',[-0.5 inf]);% shorten segment
            origTrials = 1:sEEG.trials;
            [sEEG,indx] = pop_epoch(sEEG,keyMarker,[-2 1]);% re-epoch to reponse
            rmvd = setdiff(origTrials,indx);% get removed epochs

            bEEG = pop_select(bEEG,'notrial',rmvd);% remove from baseline data
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
            cfg.padding      = 30;

            eeg = ft_freqanalysis(cfg,eeg);

            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            beeg = ft_freqanalysis(cfg,beeg);

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
                eeg.powspctrm = 10*log10(bsxfun(@rdivide,eeg.powspctrm,mean(bl.powspctrm,timedim,'omitnan'))); %  db = 10*log10(signal/baseline);
            end

            cfg = [];
            cfg.latency = [-0.4 inf];
            tmp_tfDat{cond_i,iter} = ft_selectdata(cfg,eeg);
        end
    end
    %% average iterations
    for cond_i = 1:3
        tfDat{cond_i} = ft_freqgrandaverage([],tmp_tfDat{cond_i,:});
    end
    %%
    save(outputFile,'tfDat','-v7.3')
end
