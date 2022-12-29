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
    end
end
%%
idx = tfAll.subs(:,1)==0;
tfAll.subs = tfAll.subs(~idx,:,:,:,:);
subs = subs(~idx,:);

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

%% time-resolved topopplot based on ref channel

cfg = [];
cfg.ylim = [4 10];
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
