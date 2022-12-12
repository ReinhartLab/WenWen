function singleConn(sn,IsLap,IsdePhase,IsOverWrite)
load('subs.mat');
subname = subs.name{sn};
if subs.excluded(sn)==1
    return
end
txtCell = {'','';'_lap','_dephase'};
outFile = fullfile(Dir.results,[subname,'_con',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},'.mat']);

if IsOverWrite==0 && isfile(outFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
EEG = pop_loadset('filename',set_name);
%%
EEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels
if IsLap
    X = [EEG.chanlocs.X];
    Y = [EEG.chanlocs.Y];
    Z = [EEG.chanlocs.Z];
    EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
end

EEG = pop_rmbase(EEG,[-500 -100]);% ms, remove baseline in raw

%%  re-aligned to feedback onset:29-early; 30-late; 28-correct

marker ={'S 28','S 29','S 30'};
correctCode = {{'S 28'};{'S 29','S 30'}};

for t = 1:EEG.trials
    tmp = find([EEG.event.epoch] == t);
    curEvts = {EEG.event(tmp).type};
    tmpR = find(ismember(curEvts,marker));
    if length(tmpR)>1
        EEG.event(tmp(tmpR(2:end))).type = 'rmv';
    end
end
%%
bEEG = pop_select(EEG,'time',[-1.2 1.5]); % pre-cue data

[EEG,selected]= pop_epoch(EEG,marker,[-1.2 3]);
rmvd = setdiff(1:EEG.trials,selected);

if length(rmvd)>10
    error('more than 10 trials loss after epoching')
end
if EEG.trials ~= bEEG.trials
    bEEG = pop_select(bEEG,'trial',selected);
end

%%

for cond_i = 1:2
    curCode = correctCode{cond_i};
    [sEEG,indices]=  pop_selectevent(EEG,'type',curCode);
    tmp_Trial = [EEG.event(indices).epoch];
    bsEEG = pop_select(bEEG,'trial',tmp_Trial);

    eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');
    beeg = eeglab2fieldtrip(bsEEG,'preprocessing','none');

    if IsdePhase
        avgTrl = ft_timelockanalysis([], eeg);
        eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

        avgTrl = ft_timelockanalysis([], beeg);
        beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
    end

    %%

    cfg        = [];
    cfg.method = 'wavelet';
    cfg.toi = sEEG.xmin:0.05:sEEG.xmax; % every 50ms
    cfg.output = 'fourier';
    cfg.foi = 1:40;
    cfg.width = linspace(2,7,length(cfg.foi)); % larger value for more precise f
    freq_long  = ft_freqanalysis(cfg, eeg);

    cfg.toi = bsEEG.xmin:0.05:bsEEG.xmax; % every 50ms
    freq_base = ft_freqanalysis(cfg, beeg);

    %%

    cfg           = [];
    cfg.method    = 'plv';
    plvFT = ft_connectivityanalysis(cfg, freq_long);
    plvBS = ft_connectivityanalysis(cfg, freq_base);

    cfg = [];
    cfg.latency = [-0.5 -0.1];
    tfDat{cond_i}.plvBS = ft_selectdata(cfg,plvBS);

    cfg.latency = [-0.5 0.8];
    tfDat{cond_i}.plvFT = ft_selectdata(cfg,plvFT);

end
%%
save(outFile,'tfDat','-v7.3')

