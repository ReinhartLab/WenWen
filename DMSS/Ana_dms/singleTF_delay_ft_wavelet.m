function singleTF_delay_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_delay_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

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
    clear tfDat
    ss = [1 2 4];% set size
    SStime = {[-1.6 -1.3],[-2.8 -2.5],[-5.2 -4.9]};% pretrials baseline for different SS

    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG =  pop_select(EEG,'trial',find(tmpID));

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
        cfg.pad = 20;
        eeg = ft_freqanalysis(cfg,eeg);

        if IsBL2preDelay
            cfg = [];
            cfg.baseline       = [-0.4 -0.1];% pre-delay -0.4~-0.1s as baseline
            cfg.baselinetype   = 'db';
            eeg = ft_freqbaseline(cfg, eeg);

        else %baseline to pre stim

            cfg = [];
            cfg.latency = SStime{cond_i};
            bl = ft_selectdata(cfg,eeg);

            timedim = find(size(eeg.powspctrm)==length(eeg.time));

            eeg.powspctrm = bsxfun(@(x,y)10*log10(x./y),eeg.powspctrm,mean(bl.powspctrm,timedim,'omitnan'));% dB = 10*log10(signal/baseline)
        end

        tfDat{cond_i} = eeg;

    end
    %%
    save(outputFile,'tfDat','-v7.3')
end
