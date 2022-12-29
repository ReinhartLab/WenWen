function singleTFpower(sn,IsLap,IsdePhase,IsOverWrite)
load('subs.mat');
subname = subs.name{sn};
if subs.excluded(sn)==1
    return
end
txtCell = {'','';'_lap','_dephase'};
outFile = fullfile(Dir.results,[subname,'_TFpower',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},'.mat']);

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
CodeMissing = [];
for t = 1:EEG.trials
    tmp = find([EEG.event.epoch] == t);
    curEvts = {EEG.event(tmp).type};
    tmpR = find(ismember(curEvts,marker));
    if length(tmpR)>1
        EEG.event(tmp(tmpR(2:end))).type = 'rmv'; % rename marker codes of next trial

    elseif isempty(tmpR)
        CodeMissing = [CodeMissing;t];
    end
end
EEG = pop_select(EEG,'notrial',CodeMissing);
%%
bEEG = pop_select(EEG,'time',[-1.2 1]); % pre-cue data

[EEG,selected]= pop_epoch(EEG,marker,[-1.5 3]);
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
    [sEEG,indices]=  pop_selectevent(EEG,'type',curCode);% selected EEG
    tmp_Trial = [EEG.event(indices).epoch];
    bsEEG = pop_select(bEEG,'trial',tmp_Trial);% selected baseline EEG

    eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');
    beeg = eeglab2fieldtrip(bsEEG,'preprocessing','none');

    if IsdePhase
        avgTrl = ft_timelockanalysis([], eeg);
        eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);

        avgTrl = ft_timelockanalysis([], beeg);
        beeg.trial = cellfun(@(x)x-avgTrl.avg,beeg.trial,'UniformOutput',false);
    end

    %%
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 1:40;
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.pad = 10;% in s
    cfg.toi          = beeg.time{1}(1):0.05:beeg.time{1}(end);  % the time window "slides" from -0.5 to 1.5 in 0.05 sec steps
    freq_base = ft_freqanalysis(cfg, beeg);

    %     cfg.pad = 'maxperlen';
    cfg.pad = 15;% in s
    cfg.toi          = eeg.time{1}(1):0.05:eeg.time{1}(end);  % the time window "slides" from -0.5 to 1.5 in 0.05 sec steps
    freq_long  = ft_freqanalysis(cfg, eeg);

    %%

    cfg = [];
    cfg.latency = [-0.5 1];
    tfDat{cond_i}.plvFT = ft_selectdata(cfg,freq_long);

    cfg.latency = [-0.5 -0.1];
    tfDat{cond_i}.plvBS = ft_selectdata(cfg,freq_base);

end
%%
save(outFile,'tfDat','-v7.3')

