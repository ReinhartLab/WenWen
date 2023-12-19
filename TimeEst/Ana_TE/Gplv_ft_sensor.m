addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))
clear
load('subs.mat');
txtCell = {'','','';'_lap','_dephase','_bl'};
IsLap = 0;
IsdePhase = 1;
IsBL = 0;% whether to apply baseline

for sn = 1:height(subs)
    subname = subs.name{sn};
    if subs.excluded(sn)==1
        continue
    end

    outFile = fullfile(Dir.results,[subname,'_con',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},'.mat']);
    if isfile(outFile)
        load(outFile)
    end

    if IsBL % baseline correction
        for cond_i = 1:2
            tfDat{cond_i}.plvFT.plvspctrm = bsxfun(@minus,tfDat{cond_i}.plvFT.plvspctrm,mean(tfDat{cond_i}.plvBS.plvspctrm,3));
        end
    end

    for cond_i = 1:2
        tfAll.subs(sn,cond_i,:,:,:) = tfDat{cond_i}.plvFT.plvspctrm;
        ft_subs{sn,cond_i} = tfDat{cond_i}.plvFT;
    end
end
%%
idx = tfAll.subs(:,1)==0;
tfAll.subs = tfAll.subs(~idx,:,:,:,:);
subs = subs(~idx,:);
ft_subs = ft_subs(~idx,:);

for cond_i = 1:2
    gndTF{cond_i} = tfDat{1}.plvFT;
    gndTF{cond_i}.plvspctrm = squeeze(mean(tfAll.subs(:,cond_i,:,:,:)));
end
gndTF{3} = gndTF{1};
gndTF{3}.plvspctrm = gndTF{2}.plvspctrm - gndTF{1}.plvspctrm;

condStr = {'Correct','Incorrect','Incor-Corr'};

%% time freq map of connectivity
myFigBasic
w = [1 1 1];
b = [0,0,1];
k = [0,0,0];
r = [1,0,0];
bkr = createcolormap(w,b,k,r,w); % 256x3 array

figure('Position',[100 100 1000 300]);

for cond_i = 1:3

    cfg = [];
    % cfg.xlim  = [-0.4 3];
    % cfg.ylim = [1 40];
    if cond_i ==3
        cfg.zlim = [-0.05 0.05];
    else
        %     cfg.zlim = [0 0.3];
    end
    cfg.parameter = 'plvspctrm';
    cfg.colormap = hot;
    cfg.refchannel = {'FCz'};
    cfg.channel = {'F6'};
%         cfg.channel = {'CP3'};
cfg.colormap = jet;
    cfg.figure = 'gca';

    subplot(1,3,cond_i);hold all;axis square;box on
    ft_singleplotTFR(cfg, gndTF{cond_i});

    title([condStr{cond_i}],sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]));colorbar
    ylabel('\bfFrequency(Hz)');
    xlabel('\bfTime(s)')
end
sgtitle([txtCell{IsLap+1,1}(2:end),txtCell{IsdePhase+1,2}(2:end),txtCell{IsBL+1,3}(2:end)],'FontSize',16,'Interpreter','none')

saveas(gcf,fullfile(Dir.figs,['PLV_ft_timefreqAvg_',sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.png']))
%% cluster-based permutation of incorr minus corr
cfg = [];
cfg.refchannel = {'FCz'};
cfg.channel = {'F6'};
% cfg.channel = {'CP3'};
cfg.xlim = [0.1 0.5];

tmpChancomb = ismember(tfDat{1}.plvFT.labelcmb(:,1),cfg.refchannel) & ismember(tfDat{1}.plvFT.labelcmb(:,2),cfg.channel);
tmpChancomb = find(tmpChancomb);
dat = squeeze(tfAll.subs(:,:,tmpChancomb,:,:));% sub*2*freq*time
dat = squeeze(diff(dat,1,2));%sub*freq*time
dat = permute(dat,[2 3 1]); % t*f*subs

diffstat = 'diff';%'t';% or 'diff'
nperm = 10000;
[datobs, datrnd] = cluster_test_helper(dat, nperm, diffstat);
tail = 0;alpha = .05;
clusteralpha = .05;
clusterstat = 'sum'; %'size'
[fperm.pval, fperm.h, fperm.clusterinfo] = cluster_test(datobs, datrnd, tail, alpha,...
    clusteralpha, clusterstat);           
%%
cond_i = 3;
gndTF{cond_i}.mask = repmat(fperm.h',1,1,length(gndTF{cond_i}.labelcmb));
gndTF{cond_i}.mask = permute(gndTF{cond_i}.mask ,[3 2 1]);

figure('position',[200 200 350 200]);
hold all;axis square;box on

cfg = [];
  cfg.parameter = 'plvspctrm';
  cfg.colormap = hot;
  cfg.refchannel = {'FCz'};
  cfg.channel = {'F6'};
%     cfg.channel = {'CF3'};

    cfg.colormap = jet;
  cfg.colorbartext = 'Phase locking value';
  cfg.maskparameter  = 'mask';
  cfg.maskstyle      = 'outline';
  cfg.figure = 'gca';

  ft_singleplotTFR(cfg, gndTF{cond_i});
  caxis([-0.05 0.05]);
  cb = colorbar('Ticks',[-0.05,0,0.05]);
  cb.Label.String = 'Phase locking value';
  cb.Label.Rotation = 90;

  title([condStr{cond_i}],sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]));colorbar
  ylabel('\bfFrequency(Hz)');
  xlabel('\bfTime(s)')
  saveas(gcf,fullfile(Dir.figs,['PLV_ft_timefreq_permutation_',sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.png']))
%   saveas(gcf,fullfile(Dir.figs,['PLV_ft_timefreq_permutation_',sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.pdf']))

%% averaged plv over full spectrum

myColors = crameri('bam',2);
myColors(3,:) = [0 0 0];
cfg = [];
cfg.refchannel = {'FCz'};
cfg.channel = {'F6'};
% cfg.channel = {'CP3'};
cfg.xlim = [0.1 0.5];

tmptoi = gndTF{1}.time<=cfg.xlim(2) & gndTF{1}.time>=cfg.xlim(1);
tmpChancomb = ismember(tfDat{1}.plvFT.labelcmb(:,1),cfg.refchannel) & ismember(tfDat{1}.plvFT.labelcmb(:,2),cfg.channel);
tmpChancomb = find(tmpChancomb);

subdatPLV = squeeze(mean(tfAll.subs(:,:,tmpChancomb,:,tmptoi),5));% sub*2*freq
subdatPLV(:,3,:) = diff(subdatPLV,1,2);%incor-cor

dat = squeeze(subdatPLV(:,3,:))';%freq*sub
diffstat = 'diff';%'t';% or 'diff'
nperm = 10000;

[datobs, datrnd] = cluster_test_helper(dat, nperm, diffstat);
tail = 0;alpha = .05;
clusteralpha = .05;
clusterstat = 'sum'; %'size'
[perm.pval, perm.h, perm.clusterinfo] = cluster_test(datobs, datrnd, tail, alpha,...
    clusteralpha, clusterstat);           

figure('Position',[100 100 100 200]);hold all;
rectangle('Position', [3, -1, 9, 2],'FaceColor',[ones(3,1)*0.9;0.7],'LineStyle', 'none');
for cond_i = 3%1:3
%     if cond_i == 3
%         plot(gndTF{1}.freq,squeeze(subdatPLV(:,cond_i,:)),'Color',ones(3,1)*0.8,'HandleVisibility','off')
%     end
    mn = squeeze(mean(subdatPLV(:,cond_i,:)));
    se = squeeze(std(subdatPLV(:,cond_i,:)))./sqrt(height(subs));
    plot(gndTF{1}.freq,mn,'Color',myColors(cond_i,:))
    shadedErrorBar(gndTF{1}.freq,mn,se,{'Color',myColors(cond_i,:)})

    %      hc = colorbar;
    %         hc.Label.String = sprintf('PLV (p=%.3f)',threshP);
    %         hc.Label.Rotation = -90;
    %         hc.Label.Position = [5 0 0];
    %         title(condStr{cond_i},[num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)) 'Hz: ' num2str(cfg.xlim(1)),'~',num2str(cfg.xlim(2)) 's']);
end
sig = nan(length(gndTF{1}.freq),1);
sig(perm.h) = -0.03;
plot(gndTF{1}.freq,sig,'k*','MarkerSize',2)
plot(get(gca,'xlim'),[0 0],'k:')
set(gca,'YLim',[-0.04 0.04])
xlabel('Frequency(Hz)')
ylabel('PLV');
view(90, -90); % Rotating by 90 degrees around the x-axis and -90 degrees around the y-axis
saveas(gcf,fullfile(Dir.figs,['PLV_ft_full_Spectrum',sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.png']))
saveas(gcf,fullfile(Dir.figs,['PLV_ft_full_Spectrum',sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.pdf']))
find(perm.h)' % significant freq range
%% time-resolved topopplot based on ref channel

cfg = [];
cfg.ylim = [3 12];
cfg.refchannel = 'FCz';
if IsLap
    cfg.zlim = [-0.1 0.1];
    threshP = 0.01;
else
    cfg.zlim = [-0.02 0.02];
    threshP = 0.05;
end
cfg.timeBin = {[0 0.1],[0.1 0.2],[0.2 0.3],[0.3 0.4],[0.4 0.5],[0.1 0.5]};
figure('Position',[10 10 1500 600])
for ti = 1:length(cfg.timeBin)
    cfg.xlim = cfg.timeBin{ti};

    tmpfreq = gndTF{1}.freq<=cfg.ylim(2) & gndTF{1}.freq>=cfg.ylim(1);
    tmptoi = gndTF{1}.time<=cfg.xlim(2) & gndTF{1}.time>=cfg.xlim(1);
    tmpChancomb = ismember(tfDat{1}.plvFT.labelcmb(:,1),cfg.refchannel);

    load intWM_chan_locs.mat
    chanlocs = chanlocs(ismember({chanlocs.labels},tfDat{1}.plvFT.labelcmb(tmpChancomb,2)));

    subdatTopo = tfAll.subs(:,:,tmpChancomb,tmpfreq,tmptoi);% sub*2*chan*freq*time
    subdatTopo(:,3,:,:,:) = diff(subdatTopo,1,2);%incor-cor,
    subdatTopo = squeeze(mean(mean(subdatTopo,4),5));% sub*cond*chan

    for cond_i = 1:3
        subplot(3,length(cfg.timeBin),ti+(cond_i-1)*length(cfg.timeBin));axis square

        [hmask{cond_i},pval{cond_i}]= ttest(squeeze(subdatTopo(:,cond_i,:)),0,'alpha',threshP);
%         [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval{cond_i},0.05,'pdep','yes');

        if cond_i ==3
            topoplot(squeeze(mean(subdatTopo(:,cond_i,:))),chanlocs,'plotrad',.53,'electrodes','on','style','map','emarker2',{find(hmask{cond_i}),'p','k'},'maplimits',cfg.zlim);
        else
            topoplot(squeeze(mean(subdatTopo(:,cond_i,:))),chanlocs,'plotrad',.53,'electrodes','off','style','map','emarker2',{find(hmask{cond_i}),'.','k'},'maplimits',cfg.zlim);
        end
        %                             colormap(bkr);
        hc = colorbar;
        hc.Label.String = sprintf('PLV (p=%.3f)',threshP);
        hc.Label.Rotation = -90;
        hc.Label.Position = [5 0 0];
        title(condStr{cond_i},[num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)) 'Hz: ' num2str(cfg.xlim(1)),'~',num2str(cfg.xlim(2)) 's']);
    end
end
sgtitle([txtCell{IsLap+1,1}(2:end),txtCell{IsdePhase+1,2}(2:end),txtCell{IsBL+1,3}(2:end)],'FontSize',16,'Interpreter','none')
saveas(gcf,fullfile(Dir.figs,['PLV_ft_topo_',num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)),'Hz_ref',cfg.refchannel,txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.png']))

%% plot using ft functions, cluster permutation 
cfg.refchannel = {'FCz'};
tmpChancomb = ismember(tfDat{1}.plvFT.labelcmb(:,1),cfg.refchannel);
tmpChancomb = find(tmpChancomb);
for sn = 1:length(ft_subs)
    for cond_i = 1:2
        tmp = ft_subs{sn,cond_i};
        tmp.dat = squeeze(tmp.plvspctrm(tmpChancomb,:,:));% chan*freq*time
        tmp.label = squeeze(tmp.labelcmb(tmpChancomb,2));
        tmp.dimord = 'chan_freq_time';
        ft_subs{sn,cond_i} = tmp;
    end
end

threshP = 0.05;
cfg = [];
cfg.latency          = [0.1 0.5];
cfg.frequency        = [5 12];
cfg.avgovertime = 'yes' ;
cfg.avgoverfreq = 'yes' ;
cfg.method           = 'montecarlo';
cfg.parameter = 'dat';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = threshP
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = threshP;
cfg.numrandomization = 5000;
% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb,tmp);

subj = height(subs);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, ft_subs{:,2}, ft_subs{:,1});

%% plot 
figure('Position',[100 100 300 300]);
cfg = [];
cfg.alpha  = threshP;
cfg.subplotsize = [1 1];
cfg.parameter = 'stat';
cfg.zlim   = [-5 5];
cfg.layout = 'EasycapM11.mat';
cfg.highlightsymbolseries    = ['.', '.', '+', 'o', '.'];        % the empty option will be defaulted
cfg.highlightcolorpos     = [0.8 0.8 0.8];      % the missing option will be defaulted
cfg.highlightsizeseries     = [15 15 20 20 20];      % the missing option will be defaulted
ft_clusterplot(cfg, stat);
colormap(bkr)
cb = colorbar('Ticks',[-4:4:4],'FontSize',10);
cb.Label.String = sprintf('PLV (t-value<%.3f)',threshP);
cb.Label.Rotation = -90;
cb.Label.Position = cb.Label.Position + [1 0 0];
cb.Limits = [-4 4];
title('Incor-Cor')
saveas(gcf,fullfile(Dir.figs,['PLV_ft_topo_permuted',sprintf('%.3f',threshP),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.pdf']))
