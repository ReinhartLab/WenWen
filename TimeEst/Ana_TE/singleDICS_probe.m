function singleDICS_probe(sn,IsLap,IsdePhase,IsCorretTrials)
% https://www.fieldtriptoolbox.org/tutorial/beamformer/#source_analysiswithout_contrasting_condition
load('subs.mat');
subname = subs.name{sn};

load(fullfile(Dir.cleaned,[subname,'.mat']))


sesN = 1:length(mega.EEG);
if ismember(subname,{'SE03'})
    EEG = mega.EEG{2};
else
    EEG = pop_mergeset(mega.EEG{:},sesN);
end
trials = vertcat(mega.beha{:});

clear mega
%%

EEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

cond4 = find(trials.condArray==4);
EEG = pop_select(EEG,'notrial',cond4); % remove attend trials with speed change
trials(cond4,:) = [];

if IsCorretTrials == 1
    idx = trials.rotACC~=1 | trials.rotBegin<0.25 | trials.rotBegin>3;
    trials(idx,:) = [];
    EEG = pop_select(EEG,'notrial',find(idx)); % select correct trials
end
checkEEGtrialsN(EEG,trials);

%%  re-epoch and aligned to retro-cue

ERPtime = [-4400 800];
marker = {'S201','S202'};
[EEG,indices]  = pop_epoch(EEG,marker,ERPtime./1000);
rmvEpo = setdiff(1:height(trials),indices);
trials(rmvEpo,:) = [];

if IsLap
    X = [EEG.chanlocs.X];
    Y = [EEG.chanlocs.Y];
    Z = [EEG.chanlocs.Z];
    EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
end

% change channel location into standard bem and reref to average
EEG = pop_chanedit(EEG, 'lookup',fullfile('D:\intWM-E\toolbox\eeglab2021.1\plugins\dipfit4.3\standard_BEM\elec\','standard_1020.elc'));
EEG = pop_chanedit (EEG, 'nosedir', '-Y');

EEG = pop_reref(EEG,[]); % ref to average

%%
EEG = eeglab2fieldtrip(EEG,'preprocessing','none');
if IsdePhase
    avgTrl = ft_timelockanalysis([], EEG);
    EEG.trial = cellfun(@(x)x-avgTrl.avg,EEG.trial,'UniformOutput',false);
end

% cfg = [];
% cfg.latency = [0.1 0.6];
% eeg = ft_selectdata(cfg,EEG);

cfg = [];
cfg.latency = [0 0.5];% saved as _dicsProbe05
eeg = ft_selectdata(cfg,EEG);

cfg.latency = [-0.6 -0.1]-3.8;%pre-trial baseline
eegBL = ft_selectdata(cfg,EEG);

eegcmb = ft_appenddata([],eegBL,eeg);

cfg                 = [];
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 2;
cfg.foi             = 5;
cfg.keeptrials      = 'yes';
freq_csd            = ft_freqanalysis(cfg, eegcmb);

freqBL = ft_freqanalysis(cfg,eegBL);
freqData = ft_freqanalysis(cfg,eeg);
%% read headmodel
templateheadmodel = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\headmodel\standard_bem.mat';
load(templateheadmodel); % =vol
headmodel = vol;

eeg.elec = ft_convert_units(eeg.elec, 'mm');

cfg  = [];
cfg.senstype        = 'eeg';
cfg.elec =  eeg.elec;                      % sensor positions
cfg.headmodel = headmodel;                     % volume conduction model
cfg.channel         = 'all';
cfg.grid.resolution = 1; %0.7 % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
cfg.inwardshift = 0.005; % moving sources 5 mm inwards from the skull, ...
sourcemodel = ft_prepare_leadfield(cfg);% leadfield

cfg=[];
cfg.method      = 'dics';
cfg.elec = eeg.elec;
cfg.sourcemodel = sourcemodel;       % previously computed sourcemodel
cfg.headmodel   = headmodel;        % previously computed volume conduction model
cfg.frequency   = freq_csd.freq;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '5%';% the data is rank deficient due to the low number of trials, hence some regularization is needed
cfg.dics.keepfilter  = 'yes';        % remember the filter
% cfg.dics.fixedori     = 'yes';% only keep the largest of the three dipole directions per spatial filter
cfg.dics.realfilter   = 'yes';
cfg.keepmom       = 'yes';
cfg.keepcsd       = 'yes';
sourceAll = ft_sourceanalysis(cfg, freq_csd); % contains the source estimates for all trials/both conditions

%% visualization

% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% hold all
% ft_plot_mesh(headmodel.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);%brain
% ft_plot_mesh(headmodel.bnd(2),'edgecolor','none','facealpha',0.4);%skull
% ft_plot_mesh(headmodel.bnd(3),'edgecolor','none','facealpha', 0.3,'facecolor',[0.4 0.6 0.4]);%scalp
%
% % ft_plot_sens(elec,'elec',true,'orientation',true);% electrodes
% ft_plot_sens(eeg.elec,'elec',true,'orientation',true);% electrodes
% rotate3d('on')
% view([-45 20])

%% DICS beamforming

for cond_i = 1:3

    cfg = [];
    cfg.trials = trials.condArray == cond_i;
    freqBL_Cond = ft_selectdata(cfg,freqBL);
    freqData_Cond = ft_selectdata(cfg,freqData);

    cfg=[];
    cfg.method      = 'dics';
    cfg.elec = eeg.elec;
    cfg.sourcemodel = sourcemodel;       % previously computed sourcemodel
    cfg.headmodel   = headmodel;        % previously computed volume conduction model
    cfg.frequency   = freq_csd.freq;
    % cfg.rawtrial    = 'yes';      % project each single trial through the filter. Only necessary if you are interested in reconstructing single trial data
    cfg.dics.projectnoise = 'yes';
    cfg.dics.lambda       = '5%';% the data is rank deficient due to the low number of trials, hence some regularization is needed
    cfg.dics.keepfilter  = 'yes';        % remember the filter
    % cfg.dics.fixedori     = 'yes';% only keep the largest of the three dipole directions per spatial filter
    cfg.dics.realfilter   = 'yes';
    cfg.keepmom       = 'yes';
    cfg.keepcsd       = 'yes';

    cfg.sourcemodel.filter = sourceAll.avg.filter;
    cfg.sourcemodel.label   = sourceAll.avg.label;

    source.Data{cond_i}=ft_sourceanalysis(cfg, freqData_Cond);
    source.BL{cond_i}=ft_sourceanalysis(cfg, freqBL_Cond);

end
txtCell = {'','','';'_lap','_dephase','_corrTrials'};
% save(fullfile(Dir.results,[subname,'_dicsProbe',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']),'source','-v7.3')
save(fullfile(Dir.results,[subname,'_dicsProbe05',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']),'source','-v7.3')


