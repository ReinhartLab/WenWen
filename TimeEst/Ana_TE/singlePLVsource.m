function singlePLVsource(sn,IsdePhase,IsOverWrite)
% https://www.fieldtriptoolbox.org/tutorial/networkanalysis/#connectivity-analysis-and-parcellation
% https://www.fieldtriptoolbox.org/tutorial/beamformer/#source_analysiswithout_contrasting_condition
% https://www.fieldtriptoolbox.org/tutorial/networkanalysis_eeg/
load('subs.mat');
subname = subs.name{sn};
if subs.excluded(sn)==1
    return
end
txtCell = {'','';'_lap','_dephase'};
source_conn_full.toi = [0.1 0.5];% in s, line 81&120, cfg.tapsmofrq = 2;
% source_conn_full.toi = [0.15 0.45];% in s, line 81&120, cfg.tapsmofrq  = 3;

outFile = fullfile(Dir.results,[subname,'_PLVsource',num2str(source_conn_full.toi(1)),'~',num2str(source_conn_full.toi(2)),txtCell{IsdePhase+1,2},'.mat']);

if IsOverWrite==0 && isfile(outFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
EEG = pop_loadset('filename',set_name);
%%
EEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

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
EEG= pop_epoch(EEG,marker,[-1.5 3]);

%% change channel location into standard bem and reref to average
EEG = pop_chanedit(EEG, 'lookup',fullfile('D:\intWM-E\toolbox\eeglab2021.1\plugins\dipfit4.3\standard_BEM\elec\','standard_1020.elc'));
EEG = pop_chanedit (EEG, 'nosedir', '-Y');
EEG = pop_reref(EEG,[]); % ref to average

eeg = eeglab2fieldtrip(EEG,'preprocessing','none');
if IsdePhase
    avgTrl = ft_timelockanalysis([], eeg);
    eeg.trial = cellfun(@(x)x-avgTrl.avg,eeg.trial,'UniformOutput',false);
end

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

%% get freq data

cfg = [];
cfg.latency = source_conn_full.toi;
eeg = ft_selectdata(cfg,eeg);

cfg                 = [];
cfg.output          = 'fourier';
cfg.method          = 'mtmfft';
cfg.tapsmofrq       = 2;
cfg.foi             = 6;% in Hz, given current data length & toi, exact freq varies 
freq                = ft_freqanalysis(cfg, eeg);

%% compute the beamformer filters based on the entire data.
% we call ft_sourceanalysis with ‘pcc’ as method. Essentially, this methods implements DICS (the underlying algorithm for computing the spatial filters is according to DICS),
% but provides more flexibility with respect to data handling. In this context, the advantage is that the ‘pcc’-implementation directly outputs, for each dipole location in the sourcemodel, the fourier coefficients (i.e. phase and amplitude estimates) for each of the trials.
% This can subsequently be used in a straightforward way for connectivity analysis.
% In contrast, using ‘dics’ as a method, to obtain the single trial representation of phase and amplitude is quite a bit more tedious.
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.sourcemodel       = sourcemodel;       % previously computed sourcemodel
cfg.pcc.lambda        = '10%';%how to regularise,the data is rank deficient due to the low number of trials, hence some regularization is needed
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.keeptrials        = 'yes';

sourceAll = ft_sourcedescriptives([],ft_sourceanalysis(cfg, freq));%to get the neural-activity-index

%%
curCode = correctCode{1};% correct
sEEG =  pop_selectevent(EEG,'type',curCode);% selected EEG
eegCor = eeglab2fieldtrip(sEEG,'preprocessing','none');

curCode = correctCode{2};% incorrect
sEEG =  pop_selectevent(EEG,'type',curCode);% selected EEG
eegInCor = eeglab2fieldtrip(sEEG,'preprocessing','none');

%% output fourier will give nans for plv (but only for eeg), using powandcsd instead
cfg = [];
cfg.latency = source_conn_full.toi;
eegCor = ft_selectdata(cfg,eegCor);
eegInCor = ft_selectdata(cfg,eegInCor);

cfg                 = [];
cfg.output          = 'fourier';
cfg.method          = 'mtmfft';
cfg.tapsmofrq       = 2;
cfg.foi             = 6;

freqCor             = ft_freqanalysis(cfg, eegCor);
freqInCor           = ft_freqanalysis(cfg, eegInCor);

%% use the precomputed filters
cfg                   = [];
cfg.elec = eeg.elec;
cfg.frequency         = freqCor.freq;
cfg.method            = 'pcc';
cfg.sourcemodel = sourcemodel;       % previously computed sourcemodel
cfg.headmodel         = headmodel;
cfg.pcc.lambda        = '10%';% the data is rank deficient due to the low number of trials, hence some regularization is needed
cfg.pcc.projectnoise  = 'yes';
%     cfg.pcc.keepfilter    = 'yes';
cfg.sourcemodel.filter  = sourceAll.avg.filter;

sourceCor = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freqCor)); % ft_sourcedescriptives, to get the neural-activity-index
sourceInCor = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freqInCor)); % ft_sourcedescriptives, to get the neural-activity-index

%% compute connectivity,ft_connectivityanalysis requires individual trials, and the
% combination of 'fourier' - 'pcc' is the only combination there that
% results in individual trial estimates

% compute the sparse representation,reduce memory demands and compute connectivity

sparseCor = ft_source2sparse(sourceCor);
sparseInCor = ft_source2sparse(sourceInCor);

cfg         = [];
cfg.method  ='plv';
source_connCor = ft_connectivityanalysis(cfg, sparseCor);
source_connInCor = ft_connectivityanalysis(cfg,sparseInCor);

%% then go to the 'full' representation again
source_connCor.dim     = sparseCor.dim;
% source_connCor.outside = sourceCor.outside;

source_connInCor.dim     = sparseInCor.dim;
% source_connInCor.outside = sourceInCor.outside;

source_conn_full.cor = ft_source2full(source_connCor);
source_conn_full.cor.dimord='pos_pos';

source_conn_full.incor = ft_source2full(source_connInCor);
source_conn_full.incor.dimord='pos_pos';

%%

source_conn_full.headmodel = headmodel;
source_conn_full.sourcemodel = sourcemodel;

%%
save(outFile,'source_conn_full','-v7.3')

