function singleTF_delay_ITPC_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
outputFile =fullfile(Dir.results,[subname,'_delay_ft_ITPC',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

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

        cfg = [];
        cfg.method = 'wavelet';
        cfg.foi = 1:40;
        cfg.width = linspace(2,10,length(cfg.foi)); % larger value for more precise f
        cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
        cfg.output = 'fourier';
        freq = ft_freqanalysis(cfg,eeg);

        % make a new FieldTrip-style data structure containing the ITC
        % copy the descriptive fields over from the frequency decomposition

        itc           = [];
        itc.label     = freq.label;
        itc.freq      = freq.freq;
        itc.time      = freq.time;
        itc.dimord    = 'chan_freq_time';

        F = freq.fourierspctrm;   % copy the Fourier spectrum
        N = size(F,1);           % number of trials

        % compute inter-trial phase coherence (itpc)
        itc.itpc      = F./abs(F);         % divide by amplitude
        itc.itpc      = sum(itc.itpc,1);   % sum angles
        itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
        itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension


        if IsBL2preDelay
            %             cfg = [];
            %             cfg.baseline       = [-0.4 -0.1];% pre-delay -0.4~-0.1s as baseline
            %             cfg.baselinetype   = 'db';
            %             eeg = ft_freqbaseline(cfg, eeg);

        else %baseline to pre stim

            cfg.toi = bEEG.xmin:0.05:bEEG.xmax; % every 50ms
            beeg = ft_freqanalysis(cfg,beeg);

            itcB           = [];
            itcB.label     = beeg.label;
            itcB.freq      = beeg.freq;
            itcB.time      = beeg.time;
            itcB.dimord    = 'chan_freq_time';

            F = beeg.fourierspctrm;   % copy the Fourier spectrum
            N = size(F,1);           % number of trials

            % compute inter-trial phase coherence (itpc)
            itcB.itpc      = F./abs(F);         % divide by amplitude
            itcB.itpc      = sum(itcB.itpc,1);   % sum angles
            itcB.itpc      = abs(itcB.itpc)/N;   % take the absolute value and normalize
            itcB.itpc      = squeeze(itcB.itpc); % remove the first singleton dimension

            cfg = [];
            cfg.latency = [-0.4 -0.1];
            cfg.parameter = 'itpc';
            bl = ft_selectdata(cfg,itcB);
            timedim = find(size(itcB.itpc)==length(itcB.time));

            itc.itpc = bsxfun(@(x,y)(x./y-1)*100,itc.itpc,mean(bl.itpc,timedim,'omitnan'));% (post/pre-1)*100, Arazi 2017 JN
        end

        tfDat{cond_i} = itc;

    end
    %%
    save(outputFile,'tfDat','-v7.3')

end
