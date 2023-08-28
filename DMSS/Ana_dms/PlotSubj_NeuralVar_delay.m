clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 1;
IsBL2preDelay = 1;% baselined to pre delay
IsLap = 1;
IsdePhase=1;


load('subs.mat');
txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
subname = subs.name{sn};

timeROI.all = [-5.3 4.5];% in s

SStime = {[-1.6 timeROI.all(2)],[-2.8 timeROI.all(2)],[-5.2 timeROI.all(2)]};% epoch length for different SS


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

    ss = [1 2 4];% set size
    for cond_i = 1%:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG =  pop_select(EEG,'trial',find(tmpID));

        bEEG = pop_select(sEEG,'time',[-inf 2]);% re-epoch into shorter segments
        timelimits = [-0.88 1.5];% pre-stimuli baseline
        [bEEG,indices] = pop_epoch(bEEG,stim1,timelimits);

        rmv = setdiff(1:sEEG.trials,indices);
        sEEG = pop_select(sEEG,'notrial',rmv);

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

        if IsBL2preDelay==0 %baseline to pre trial
            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            beeg = ft_freqanalysis(cfg,beeg);
        end

        if IsBL2preDelay
            cfg = [];
            cfg.latency = [-0.4 -0.1];% select baseline as preDelay
            bl = ft_selectdata(cfg,eeg);

        else %baseline to pre trial

            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,beeg);% select baseline as pre trial
        end

        %%
        freqID.band = [15 25];
        freqID.idx = dsearchn(eeg.freq',freqID.band');
        freqID.idx = freqID.idx(1):freqID.idx(2);

        chanROI = {'Fz','F1','F2','FCz','AFz'};
        clear chanID
        [~,chanID] = ismember(chanROI,eeg.label);
        %%
        clear dat
        timedim = find(size(eeg.powspctrm)==length(eeg.time));
        dat.nv = squeeze(log(var(eeg.powspctrm,0,1,"omitnan")./mean(var(bl.powspctrm,0,1,"omitnan"),timedim)));%chan*freq*time

        av = ft_freqdescriptives([],eeg);% avg

        cfg = [];
        cfg.baselinetype = 'db';
        cfg.baseline = [-0.4 -0.1];
        av = ft_freqbaseline(cfg,av);
        eeg = ft_freqbaseline(cfg,eeg);

        dat.trials = squeeze(mean(mean(eeg.powspctrm(:,chanID,freqID.idx,:),3),2));% trial*time

        %%
        myFigBasic
        figure('Position',[100 100 250 250]);hold all
        plot(eeg.time,dat.trials','Color',ones(3,1)*0.9,'HandleVisibility','off')
        plot(eeg.time,squeeze(mean(mean(dat.nv(chanID,freqID.idx,:),2),1)),'LineWidth',2)
        plot(eeg.time,squeeze(mean(mean(av.powspctrm(chanID,freqID.idx,:),2),1)),'LineWidth',2)
        legend({'Mean','Relative variance'})
        ylabel('Amplitude in db or a.u.')
        xlabel('Time (0s=Maintenance)')
        plot([0 0],get(gca,'YLim'),'k--','LineWidth',0.1,'HandleVisibility','off')
        plot([3 3],get(gca,'YLim'),'k--','LineWidth',0.1,'Color',[1 1 1]*0.1,'HandleVisibility','off')
        set(gca,'XTick',[-4.8 -3.6 -2.4 -1.2 0 3 timeROI.all(2)],'XLim',SStime{cond_i})
        title(subname)
        saveas(gcf,fullfile(Dir.figs,'IndividualResults',[subname,'NeuVar_delay_load',num2str(cond_i),[chanROI{:}],txtCell{IsLap+1,1},num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
        saveas(gcf,fullfile(Dir.figs,'IndividualResults',[subname,'NeuVar_delay_load',num2str(cond_i),[chanROI{:}],txtCell{IsLap+1,1},num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))
    end
end
