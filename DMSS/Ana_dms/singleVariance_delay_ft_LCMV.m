function singleVariance_delay_ft_LCMV(sn,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_delay_LCMV',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

if IsOverwrite==0 && isfile(outputFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
if isfile(set_name)
    EEG = pop_loadset('filename',set_name);% baseline corrected to pre-stim
    EEG = pop_select(EEG,'nochannel',{'LHEOG','RHEOG','BVEOG','TVEOG'});
      
    %% baseline correction: For a correct covariance estimation it is important that you used the cfg.demean = yes when the function ft_preprocessing was applied.

    if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
        stim1 = {'S 11','S 21', 'S 31'};
    else
        stim1 = {'11','21','31'};
    end

    bEEG = pop_select(EEG,'time',[-inf 2]);% re-epoch into shorter segments
    timelimits = [-0.88 1];% pre-stimuli baseline
    [bEEG,indices]= pop_epoch(bEEG,stim1,timelimits);
    rmv = setdiff(1:EEG.trials,indices);
    EEG = pop_select(EEG,'notrial',rmv);
   
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
 %% change channel location into standard bem and reref to average: https://www.fieldtriptoolbox.org/tutorial/beamformer_lcmv/
    EEG = pop_chanedit(EEG, 'lookup',fullfile('D:\intWM-E\toolbox\eeglab2021.1\plugins\dipfit4.3\standard_BEM\elec\','standard_1020.elc'));
    EEG = pop_chanedit (EEG, 'nosedir', '-Y');
    EEG = pop_reref(EEG,[]); % ref to average

    %%
    if IsBL2preDelay
        EEG = pop_rmbase(EEG,[-200 0]);% in ms
    else
        tmp_toi = bEEG.times<=0 & bEEG.times>=-200;
        EEG.data = bsxfun(@minus,EEG.data,mean(bEEG.data(:,tmp_toi,:),2));
    end

    eeg = eeglab2fieldtrip(EEG,'preprocessing','none');
    beeg = eeglab2fieldtrip(bEEG,'preprocessing','none');

    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [15 25];% in Hz,https://www.fieldtriptoolbox.org/tutorial/virtual_sensors/
    eeg = ft_preprocessing(cfg,eeg);
    beeg =  ft_preprocessing(cfg,beeg);

    cfg                  = [];
    cfg.covariance       = 'yes';
    if IsCorretTrials==1
        cfg.trials = M.button_resp_corr ==1;
    end
    timelock             = ft_timelockanalysis(cfg, eeg);
    timelock_b           = ft_timelockanalysis(cfg, beeg);

    %% read headmodel
    templateheadmodel = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\headmodel\standard_bem.mat';
    load(templateheadmodel); % =vol
    headmodel = vol;

    eeg.elec = ft_convert_units(eeg.elec, 'mm');

    cfg  = [];
    cfg.elec =  eeg.elec;                      % sensor positions
    cfg.headmodel = headmodel;                     % volume conduction model
    cfg.grid.resolution = 1; % use a 3-D grid with a 1 cm resolution
    cfg.grid.unit       = 'cm';
    cfg.inwardshift = 0.005; % moving sources 5 mm inwards from the skull, ...
    sourcemodel = ft_prepare_leadfield(cfg);%

    %% compute the beamformer filters based on the entire data.

    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.sourcemodel      = sourcemodel; % leadfield
    cfg.headmodel        = headmodel; % volume conduction model (headmodel)
    cfg.lcmv.keepfilter  = 'yes';
    cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
    sourceAll            = ft_sourceanalysis(cfg, timelock);
    sourceAll_b          = ft_sourceanalysis(cfg, timelock_b);

    %%
    ss = [1 2 4];% set size
    for cond_i = 1:3
        tmpID = M.ss_num == ss(cond_i);
        if IsCorretTrials ==1
            tmpID = M.ss_num == ss(cond_i) & M.button_resp_corr ==1;
        end
        sEEG =  pop_select(EEG,'trial',find(tmpID));
        seeg = eeglab2fieldtrip(sEEG,'preprocessing','none');

        blEEG =  pop_select(bEEG,'trial',find(tmpID));
        bleeg = eeglab2fieldtrip(blEEG,'preprocessing','none');

        cfg                  = [];
        cfg.covariance       = 'yes';
        timelock_s           = ft_timelockanalysis(cfg, seeg);
        timelock_bl          = ft_timelockanalysis(cfg, bleeg);

        cfg                  = [];
        cfg.method           = 'lcmv';
        cfg.sourcemodel      = sourcemodel; % leadfield
        cfg.headmodel        = headmodel; % volume conduction model (headmodel)
        cfg.lcmv.keepfilter  = 'yes';
        cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
        cfg.projectnoise     = 'yes';
        %  neural activity index (NAI), remove the center of the head bias shown above.
        %         The NAI is the power normalized with an estimate of the spatially inhomogeneous noise.
        % by setting cfg.dics.projectnoise=’yes’
        cfg.sourcemodel.filter  = sourceAll.avg.filter;

        source{cond_i}          = ft_sourceanalysis(cfg, timelock_s);

      

        if IsBL2preDelay==0
            cfg.sourcemodel.filter  = sourceAll_b.avg.filter;
            source_b{cond_i}        = ft_sourceanalysis(cfg, timelock_bl);

            beamformerBL = vertcat(source_b{cond_i}.avg.filter{:});%n dipole * n chan
            source_BL = bleeg;
            source_BL = rmfield(source_BL,{'cfg','trialinfo','label','elec'});
            source_BL.label = cellstr(num2str(find(sourceAll_b.inside)));

            for i=1:length(bleeg.trial)
                source_BL.trial{i} = beamformerBL * bleeg.trial{i};
            end
        end

        beamformer = vertcat(source{cond_i}.avg.filter{:});%n dipole * n chan
        source_data = seeg;
        source_data = rmfield(source_data,{'cfg','trialinfo','label','elec'});
        source_data.label = cellstr(num2str(find(sourceAll.inside)));

        for i=1:length(seeg.trial)
            source_data.trial{i} = beamformer * seeg.trial{i};
        end

        if IsdePhase
            avgTrl = ft_timelockanalysis([], source_data);
            source_data.trial = cellfun(@(x)x-avgTrl.avg,source_data.trial,'UniformOutput',false);

            if IsBL2preDelay==0
                avgTrl = ft_timelockanalysis([], source_BL);
                source_BL.trial = cellfun(@(x)x-avgTrl.avg,source_BL.trial,'UniformOutput',false);
            end
        end

        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.out = 'pow';
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.keeptrials  = 'yes';
        cfg.pad = 30;
        source_freq = ft_freqanalysis(cfg,source_data);

        if IsBL2preDelay==0 %baseline to pre trial
            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            source_freqBL = ft_freqanalysis(cfg,source_BL);
        end

        if IsBL2preDelay
            cfg = [];
            cfg.latency = [-0.4 -0.1];% select baseline as preDelay
            bl = ft_selectdata(cfg,source_freq);

        else %baseline to pre trial
            cfg = [];
            cfg.latency = [-0.4 -0.1];
            bl = ft_selectdata(cfg,source_freqBL);% select baseline as pre trial
        end

        % variance
        timedim = find(size(source_freq.powspctrm)==length(source_freq.time));
        nv = squeeze(log(var(source_freq.powspctrm,0,1,'omitnan')./mean(var(bl.powspctrm,0,1,'omitnan'),timedim)));%chan*freq*time

        sourceVar = ft_freqdescriptives([],source_freq);% get ft structure;
        sourceVar.powspctrm = nv;

        tfDat{cond_i} = sourceVar;
    end

    %%
    save(outputFile,'tfDat','sourcemodel','headmodel','-v7.3')
end
