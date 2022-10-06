function [sn,threshV,locthresh,globthresh,isPlot] = genMega(sn,varargin)

p=inputParser;

default_threshV = 75;
default_locthresh = 8;
default_globthresh = 8;
default_zthresh = 8;

default_isPlot = 1;

addParameter(p,'threshV',default_threshV,@isnumeric);
addParameter(p,'locthresh',default_locthresh,@isnumeric);
addParameter(p,'globthresh',default_globthresh,@isnumeric);
addParameter(p,'isPlot',default_isPlot,@isnumeric);
addParameter(p,'zthresh',default_zthresh,@isnumeric);


parse(p,varargin{:});

threshV = p.Results.threshV;
locthresh = p.Results.locthresh;
globthresh = p.Results.globthresh;
isPlot = p.Results.isPlot;
zthresh = p.Results.zthresh;
% isPlot =1;
% threshV = 100;
% locthresh = 8;
% globthresh = 8;

%%

load('subs.mat');
subname = subs.name{sn};

set_name{1} = fullfile(Dir.prepro,[subname,'_postICA.set']);
behaFile{1} = fullfile(Dir.beha,subname,[subname '.mat']);

ses2 = fullfile(Dir.raw,[subname,'02']);
if isfolder(ses2)
    sesN = 2;
    set_name{2} = fullfile(Dir.prepro,[subname,'02','_postICA.set']);
    behaFile{2} = fullfile(Dir.beha,[subname '02'],[subname '02.mat']);
else
    sesN = 1;
end

matname = 'trialRM.mat';
if exist(matname)
    load(matname)
end

%%
for ses_i = 1:sesN
    EEG = pop_loadset('filename',set_name{ses_i});

    load(behaFile{ses_i});
    if ismember('pause',trials.Properties.VariableNames)
        trials.pause = [];% remove column, some ppts don't have 'pause'
    end

    if EEG.trials ~= height(trials)
        if isfield(trialRM,'Epoch')
            trials(trialRM.Epoch{ses_i,sn},:) = [];
        end% trials removed due to discontinuity
        if isfield(trialRM,'artifact')
            trials(trialRM.artifact{ses_i,sn},:) = [];% trials manually removed by visual inspection
        end
    end

    checkEEGtrialsN(EEG,trials);

    %%
    sEEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'},'time',[-0.5 4.5]);% remove EOG channels

    [~, rejT] = pop_eegthresh( sEEG, 1, 1:sEEG.nbchan, -threshV, threshV, sEEG.xmin, sEEG.xmax,0,0);
    [~, ~, ~, rejK] = pop_rejkurt_ww( sEEG, 1, 1:sEEG.nbchan, locthresh,globthresh, 'reject',0,'plotflag',0);
    [~, ~, ~, rejP] = pop_jointprob_ww( sEEG, 1, 1:sEEG.nbchan, locthresh, globthresh,1, 'reject',0,'plotflag',0);

    tmp = abs(zscore(sEEG.data,0,3))>zthresh;% along trial dimension
    rejZ = squeeze(sum(sum(tmp),2))>0;

    tmp = abs(zscore(sEEG.data,0,2))>zthresh;% along time dimension
    rejZ2 = squeeze(sum(sum(tmp),2))>0;
    %%
    rej_union = union(rejT,[find(rejK) find(rejP) find(rejZ)' find(rejZ2)']);
    fprintf('\rRejecting %d trials\n',length(rej_union))

    if isPlot
        rEEG = pop_select(EEG,'trial',rej_union);
        pop_eegplot(rEEG,1,1,1);
    end

    trialRM.postRejction{ses_i,sn} = rej_union;
    save(matname,'trialRM')

    EEG = pop_rejepoch( EEG, rej_union,0);
    trials(rej_union,:)= [];
    checkEEGtrialsN(EEG,trials);

    if isPlot
        rEEG = pop_select(sEEG,'notrial',rej_union);
        pop_eegplot(rEEG,1,1,1);
    end

    mega.EEG{ses_i} = EEG;
    mega.beha{ses_i} = trials;

end
fprintf('Saving mega file for %s\r',subname)
megaName = fullfile(Dir.cleaned,[subname,'.mat']);
save(megaName,'mega','-v7.3')
