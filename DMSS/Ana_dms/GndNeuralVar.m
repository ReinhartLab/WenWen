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
IsBL2preDelay = 0;
IsLap = 0;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Young-Old'};
condStr = {'Load 1','Load 2','Load 4'};
condDiffStr = {'S2-S1','S4-S1','S4-S2'};

% frontalROI = {'Cz','C1','C2','CPz'};
% frontalROI = {'Fz','F1','F2'};
frontalROI = {'Fz','F1','F2','FCz','FC1','FC2'};
occipROI = {'Pz','POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz
timeROI = [-1 4.5];% in s
SStime = {[-1.6 4.5],[-2.8 4.5],[-5.2 4.5]};% epoch length for different SS

timeROIdelayFtest = [-1 3];% in s
timeROIprobeFtest = [3.0 3.6];% in s

%%
subsAll = cell(subN,3);
for sub_i = 1:subN
    subname = subs.name{sub_i};
    matName = fullfile(Dir.results,[subname,'_var_delay_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);
        subsAll(sub_i,:) = tfDat;
    end
end

empty_idx = cellfun(@isempty,subsAll(:,1));
subsAll(empty_idx,:) = [];
subs(empty_idx,:) = [];

for gi = 1:2
    snList = subs.group == gi;
    for ss = 1:3
        gndTF{gi,ss} = ft_freqgrandaverage([],subsAll{snList,ss});

        cfg = [];
        cfg.keepindividual = 'yes';
        tmp = ft_freqgrandaverage(cfg,subsAll{snList,ss});
        gndTF{gi,ss}.indv = tmp.powspctrm;
    end
    GroupAvg{gi} = ft_freqgrandaverage([],subsAll{snList,:});
end

for gi = 1:2
    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'powspctrm';
    gndTF_SSdiff{gi,1} = ft_math(cfg, gndTF{gi,2},gndTF{gi,1}); %set size difference
    gndTF_SSdiff{gi,2} = ft_math(cfg, gndTF{gi,3},gndTF{gi,1});
    gndTF_SSdiff{gi,3} = ft_math(cfg, gndTF{gi,3},gndTF{gi,2});
end

for ss = 1:3
    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'powspctrm';
    gndTF_Groupdiff{ss} = ft_math(cfg, gndTF{1,ss},gndTF{2,ss}); %group difference
end

GroupAvg{3}  = ft_math(cfg, GroupAvg{1},GroupAvg{2});
%% time-frequency map
myFigBasic;

w = [1 1 1];
b = [0,0,1];
k = [0,0,0];
r = [1,0,0];
bkr = createcolormap(w,b,k,r,w); % 256x3 array

for gi = 1:2
    figure('Position',[100 100 1000 400]);

    for cond_i = 1:3
        subplot(2,3,cond_i);
        cfg = [];
        cfg.xlim  = SStime{cond_i};
        cfg.ylim = [1 40];
        cfg.zlim = [-1 1];
        % cfg.colormap = bkr;
        cfg.colormap = jet;

        % cfg.channel = {'all'};
        cfg.channel = frontalROI;
        cfg.figure = 'gca';
        % cfg.maskparameter = 'mask';
        % cfg.maskstyle = 'outline';

        % cfg.maskstyle = 'opacity';
        % cfg.maskalpha = 0.3;

        ft_singleplotTFR(cfg, gndTF{gi,cond_i});

        title([condStr{cond_i},' ' groupStr{gi}],'FontSize',14);
        hc = colorbar;
        hc.Label.String = 'Relative variance(AU)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',[-4.8 -3.6 -2.4 -1.2 0 3])
        if cond_i == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time(0s=maintenance)')
        rectangle('Position',[0 0.5 3 40],'LineWidth',1,'LineStyle',':')
        rectangle('Position',[-1 0.5 4 40],'LineWidth',1,'LineStyle',':')

        subplot(2,3,3+cond_i);
        cfg = [];
        % cfg.ylim = [15 25];
        cfg.ylim = freq.betaFreq;%freq

        cfg.xlim = timeROI;%time
        cfg.zlim = [-0.6 0];
        cfg.highlightchannel = frontalROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';

        ft_topoplotTFR(cfg, gndTF{gi,cond_i});colorbar;

    end
    saveas(gcf,fullfile(Dir.figs,['NeuVar_',groupStr{gi},[frontalROI{:}],txtCell{IsLap+1,1},num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))
end

%%  time frequency of group average
cfg = [];
cfg.xlim  = [-1 5];
cfg.ylim = [1 40];
cfg.zlim = [-1 1];
% cfg.colormap = bkr;
cfg.colormap = jet;
cfg.channel = frontalROI;
cfg.figure = 'gca';
% cfg.maskparameter = 'mask';
% cfg.maskstyle = 'outline';
% cfg.maskstyle = 'opacity';
% cfg.maskalpha = 0.3;

figure('Position',[100 100 800 290]);
for gi = 1:3
    subplot(1,3,gi);axis square;hold all

    ft_singleplotTFR(cfg, GroupAvg{gi});

    title(groupStr{gi},'Average set sizes','FontSize',14);
    hc = colorbar;
    hc.Label.String = 'Relative variance(AU)';
    hc.Label.Rotation = -90;
    hc.Label.Position = [4.8 0 0];

    xlabel('\bfTime(s)')
    rectangle('Position',[0 0.5 3 40],'LineWidth',1,'LineStyle','--')
    set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',[-1 0 3],'YLim',[0.5 40])
end
saveas(gcf,fullfile(Dir.figs,['NeuVar_GroupAvg_',[cfg.channel{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))
%% time frequency of group diff

figure('Position',[100 100 800 400])
for cond_i = 1:3
    subplot(2,3,cond_i);axis square
    cfg = [];
    cfg.xlim  = [-1 4.5];
    cfg.ylim = [1 40];
    cfg.zlim = [-1 1];
    % cfg.colormap = bkr;
    cfg.colormap = jet;

    cfg.channel = frontalROI;
    cfg.figure = 'gca';
    % cfg.maskparameter = 'mask';
    % cfg.maskstyle = 'outline';
    % cfg.maskstyle = 'opacity';
    % cfg.maskalpha = 0.3;

    ft_singleplotTFR(cfg, gndTF_Groupdiff{cond_i});
    set(gca,'xtick',[-1 0 3])
    rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')
    title(condStr{cond_i},groupStr{3},'FontSize',14);colorbar

    subplot(2,3,cond_i+3);
    cfg = [];
    cfg.ylim = freq.betaFreq;%freq
    cfg.xlim = timeROIdelayFtest;%time
    cfg.zlim = [-0.8 0.8];
    cfg.highlightchannel = frontalROI;
    cfg.layout = 'easyCapM1';
    cfg.colormap = jet;
    % cfg.marker = 'labels';
    cfg.markersymbol = '.';
    cfg.highlight = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 18;
    cfg.figure      = 'gca';
    ft_topoplotTFR(cfg, gndTF_Groupdiff{cond_i});
    title(condStr{cond_i},groupStr{3});colorbar

end
saveas(gcf,fullfile(Dir.figs,['NeuVarGroupDiff_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))

%%
freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID = dsearchn(gndTF{1}.time',timeROI');
timeID = timeID(1):timeID(2);

timeIDdelayFtest = dsearchn(gndTF{1}.time',timeROIdelayFtest');
timeIDdelayFtest = timeIDdelayFtest(1):timeIDdelayFtest(2);

timeIDprobeFtest = dsearchn(gndTF{1}.time',timeROIprobeFtest');
timeIDprobeFtest = timeIDprobeFtest(1):timeIDprobeFtest(2);

clear chanID
[log1,chanID.frontal] = ismember(frontalROI,gndTF{1}.label);
[log2,chanID.occip] = ismember(occipROI,gndTF{1}.label);

clear dat SStimeID

for gi = 1:2
    for cond_i = 1:3
        SStimeID{cond_i} = dsearchn(gndTF{1}.time',SStime{cond_i}');
        SStimeID{cond_i} = SStimeID{cond_i}(1):SStimeID{cond_i}(2);

        dat.betaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeIDdelayFtest),2),3),4));
        dat.betaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeIDdelayFtest),2),3),4));
        dat.betaAvgFrontalProbe{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeIDprobeFtest),2),3),4));

        dat.alphaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeIDdelayFtest),2),3),4));
        dat.alphaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeIDdelayFtest),2),3),4));

        dat.betaCurveOccip{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,SStimeID{cond_i}),2),3));
        dat.betaCurveFrontal{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,SStimeID{cond_i}),2),3));

        dat.alphaCurveOccip{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,SStimeID{cond_i}),2),3));
        dat.alphaCurveFrontal{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,SStimeID{cond_i}),2),3));
    end
end

% mixed ANOVA: group*cond = 2*3
clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontal{1},size(dat.betaAvgFrontal{1},1)*size(dat.betaAvgFrontal{1},2),1);...
    reshape(dat.betaAvgFrontal{2},size(dat.betaAvgFrontal{2},1)*size(dat.betaAvgFrontal{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.betaFrontal, DFs.betaFrontal, MSQs.betaFrontal, Fs.betaFrontal, Ps.betaFrontal]=mixed_between_within_anova(Xa);

clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontalProbe{1},size(dat.betaAvgFrontalProbe{1},1)*size(dat.betaAvgFrontalProbe{1},2),1);...
    reshape(dat.betaAvgFrontalProbe{2},size(dat.betaAvgFrontalProbe{2},1)*size(dat.betaAvgFrontalProbe{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';
[SSQs.betaFrontalProbe, DFs.betaFrontalProbe, MSQs.betaFrontalProbe, Fs.betaFrontalProbe, Ps.betaFrontalProbe]=mixed_between_within_anova(Xa);

clear Xa
Xa(:,1) = [reshape(dat.betaAvgOccip{1},size(dat.betaAvgOccip{1},1)*size(dat.betaAvgOccip{1},2),1);...
    reshape(dat.betaAvgOccip{2},size(dat.betaAvgOccip{2},1)*size(dat.betaAvgOccip{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.betaOccip, DFs.betaOccip, MSQs.betaOccip, Fs.betaOccip, Ps.betaOccip]=mixed_between_within_anova(Xa);

clear Xa
Xa(:,1) = [reshape(dat.alphaAvgOccip{1},size(dat.alphaAvgOccip{1},1)*size(dat.alphaAvgOccip{1},2),1);...
    reshape(dat.alphaAvgOccip{2},size(dat.alphaAvgOccip{2},1)*size(dat.alphaAvgOccip{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.alpha, DFs.alpha, MSQs.alpha, Fs.alpha, Ps.alpha]=mixed_between_within_anova(Xa);

%% linear fit of SS
clear line_fit line_fitProbe
for gi = 1:2
    for sub_i = 1:size(dat.betaAvgFrontal{gi},1)
        line_fit{gi}.params(sub_i,:) = polyfit([1 2 4]',dat.betaAvgFrontal{gi}(sub_i,:)',1);
        line_fitProbe{gi}.params(sub_i,:) = polyfit([1 2 4]',dat.betaAvgFrontalProbe{gi}(sub_i,:)',1);
    end
    [~,line_fit{gi}.p,~,line_fit{gi}.Tstat] = ttest(line_fit{gi}.params(:,1));
    [~,line_fitProbe{gi}.p,~,line_fitProbe{gi}.Tstat] = ttest(line_fitProbe{gi}.params(:,1));

end

%% bar & time series
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
figure('Position',[100 100 800 1200]);

subplot(4,2,5);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
plot(xcord(1,:),dat.betaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

for gi = 1:2
    if line_fit{gi}.p<.05
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'r','LineWidth',1.5,'HandleVisibility','off')
    else
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
    end
end
legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')
title(sprintf('%s beta\ndelay (%.1f~%.1fs)',[frontalROI{:}],timeROIdelayFtest(1),timeROIdelayFtest(2)))
% plot significance
if Ps.betaFrontal{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontal{gi} dat.betaAvgFrontal{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-0.8,-0.6,3);
        xposi = {[1 2],[2 3],[1 3]};
        for ss = 1:3
            if tmp_pval(ss)<=.05
                plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
                text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
            end
        end
    end
end

subplot(4,2,6);hold all;axis square
mn = cellfun(@mean,dat.betaAvgOccip,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgOccip,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.betaAvgOccip{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgOccip{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')
title(sprintf('Posterior beta\n delay (%.1f~%.1fs)',timeROIdelayFtest(1),timeROIdelayFtest(2)))

subplot(4,2,7);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontalProbe,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontalProbe,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.betaAvgFrontalProbe{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontalProbe{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';

for gi = 1:2
    if line_fitProbe{gi}.p<.05
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'r','LineWidth',1.5,'HandleVisibility','off')
    else
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
    end
end

legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')
title(sprintf('%s beta\nprobe (%.1f~%.1fs)',[frontalROI{:}],timeROIprobeFtest(1),timeROIprobeFtest(2)))
if Ps.betaFrontal{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontalProbe{gi} dat.betaAvgFrontalProbe{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-1.3,-1,3);
        xposi = {[1 2],[2 3],[1 3]};
        for ss = 1:3
            if tmp_pval(ss)<=.05
                plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
                text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
            end
        end
    end
end

myColors = crameri('batlowW',5);

% times = gndTF{1}.time(timeID);
for gi = 1:2
    subplot(4,2,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveFrontal{gi}{cond_i}));
        se = squeeze(std(dat.betaCurveFrontal{gi}{cond_i},0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(gndTF{1}.time(SStimeID{cond_i}),mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([-4.8 -4.8],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-3.6 -3.6],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-2.4 -2.4],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-1.2 -1.2],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    plot([3 3],get(gca,'YLim'),'k--','HandleVisibility','off')

    legend(condStr)
    ytickformat('%.1f')
    ylabel('Variance(AU)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[frontalROI{:}],freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XTick',[-4.8 -3.6 -2.4 -1.2 0 3 4.5])

    subplot(4,2,gi+2);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveOccip{gi}{cond_i}));
        se = squeeze(std(dat.betaCurveOccip{gi}{cond_i},0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(gndTF{1}.time(SStimeID{cond_i}),mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    plot([3 3],get(gca,'YLim'),'k--','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Variance(AU)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI,'XTick',[-4.8 -2.4 -1.2 0 3 4.5])
end

saveas(gca,fullfile(Dir.figs,['DelayBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))
%%
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[100 100 800 1000]);

subplot(3,2,5);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgOccip,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgOccip,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
plot(xcord(1,:),dat.alphaAvgOccip{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgOccip{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')
title(sprintf('Posterior alpha(%.1f~%.1fs)',timeROIdelayFtest(1),timeROIdelayFtest(2)))

subplot(3,2,6);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.alphaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')

title(sprintf('%s alpha(%.1f~%.1fs)',[frontalROI{:}],timeROIdelayFtest(1),timeROIdelayFtest(2)))

myColors = crameri('batlowW',5);

times = gndTF{1}.time(timeID);
for gi = 1:2
    subplot(3,2,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveOccip{gi}{cond_i}));
        se = squeeze(std(dat.alphaCurveOccip{gi}{cond_i},0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(gndTF{1}.time(SStimeID{cond_i}),mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([-4.8 -4.8],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-3.6 -3.6],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-2.4 -2.4],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-1.2 -1.2],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Variance(AU)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.alphaFreq(1),freq.alphaFreq(2)))

    subplot(3,2,gi+2);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveFrontal{gi}{cond_i}));
        se = squeeze(std(dat.alphaCurveFrontal{gi}{cond_i},0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(gndTF{1}.time(SStimeID{cond_i}),mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([-4.8 -4.8],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-3.6 -3.6],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-2.4 -2.4],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-1.2 -1.2],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Variance(AU)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[frontalROI{:}],freq.alphaFreq(1),freq.alphaFreq(2)))
end

saveas(gca,fullfile(Dir.figs,['DelayAlphaVariance_',num2str(freq.alphaFreq(1)),'~',num2str(freq.alphaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))

%% GLM beta coeff
clear glmBeta
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);
dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));

myColors = crameri('bamako',4);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[500 500 650 600]);
threshP = .05;
for gi = 1:2
    subplot(2,2,gi);hold all;axis square;
    for cond_i = 1:3
        X = dat.beha{gi}(:,cond_i);
        for t = 1:size(dat.betaCurveFrontal{gi}{cond_i},2)
            Y = squeeze(dat.betaCurveFrontal{gi}{cond_i}(:,t));
            [b,dev,stats] = glmfit(X,Y,'normal');
            glmBeta.b{gi,cond_i}(t) = b(2);
            glmBeta.p{gi,cond_i}(t) = stats.p(2);
        end

        plot(gndTF{1}.time(SStimeID{cond_i}),glmBeta.b{gi,cond_i},'Color',myColors(cond_i,:))

        mySig = nan(length(gndTF{1}.time(SStimeID{cond_i})),1);
        mySig(glmBeta.p{gi,cond_i}(t)<threshP) = min(glmBeta.b{gi,cond_i})*1.1;
        plot(gndTF{1}.time(SStimeID{cond_i}),mySig,'Color',myColors(cond_i,:),'LineWidth',3,'HandleVisibility','off')
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('\beta coefficient(variance~accuracy)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},'Frontal beta')

    subplot(2,2,gi+2);hold all;axis square;

    for cond_i = 1:3
        X = dat.beha{gi}(:,cond_i);
        for t = 1:size(dat.betaCurveOccip{gi}{cond_i},2)
            Y = squeeze(dat.betaCurveOccip{gi}{cond_i}(:,t));
            [b,dev,stats] = glmfit(X,Y,'normal');
            glmBeta.b{gi,cond_i}(t) = b(2);
            glmBeta.p{gi,cond_i}(t) = stats.p(2);
        end

        plot(gndTF{1}.time(SStimeID{cond_i}),glmBeta.b{gi,cond_i},'Color',myColors(cond_i,:))

        mySig = nan(length(gndTF{1}.time(SStimeID{cond_i})),1);
        mySig(glmBeta.p{gi,cond_i}(t)<threshP) = min(glmBeta.b{gi,cond_i})*1.1;
        plot(gndTF{1}.time(SStimeID{cond_i}),mySig,'Color',myColors(cond_i,:),'LineWidth',3,'HandleVisibility','off')
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('\beta coefficient(variance~accuracy)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},'Posterior beta')

end
saveas(gca,fullfile(Dir.figs,['DelayBetaVarBehavGLM_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);

behaStr = {'RT','acc'};

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',5);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

for idx = 1:2
    if idx == 2
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));
    else
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_rt','s2_rt','s3_rt'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_rt','s2_rt','s3_rt'}));
    end

    clear g
    for gi = 1:2
        [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Spearman');
        %     [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Pearson');

        for cond = 1:3
            Y = dat.beha{gi}(:,cond);
            X = dat.betaAvgFrontal{gi}(:,cond);

            g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
            g(gi,cond).geom_point(); %Draw a scatter plot
            g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
            g(gi,cond).set_names('x','Beta variance','y',behaStr{idx}); %Set the title of the axis
            g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
            g(gi,cond).set_title(sprintf('%s: Rho = %.3f, p = %.3f',groupStr{gi},r(cond,cond),pval(cond,cond)));
            g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
        end
        %%

        clear X Y
        Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
        X(:,1) = sort(repmat([1 2 4]',size(dat.betaAvgFrontal{gi},1),1));
        X(:,2) = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);
        mdl_124 =  stepwiseglm(X,Y,'linear','Distribution','normal','Criterion','aic','Lower','linear');

        clear X Y
        Y = reshape(dat.beha{gi}(:,[2 3]),size(dat.beha{gi},1)*2,1);
        X(:,1) = sort(repmat([2 4]',size(dat.betaAvgFrontal{gi},1),1));
        X(:,2) = reshape(dat.betaAvgFrontal{gi}(:,[2 3]),size(dat.betaAvgFrontal{gi},1)*2,1);
        mdl_24 =  stepwiseglm(X,Y,'linear','Distribution','normal','Criterion','aic','Lower','linear');

        %%
        cond = 4;
        Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
        X = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);
        condArray = sort(repmat([1 2 3]',sum(subs.group==gi),1));
        g(gi,cond) = gramm('x',X,'y',Y,'color',condStr(condArray)); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,cond).geom_point(); %Draw a scatter plot
        g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
        g(gi,cond).set_names('x','Beta variance','y',behaStr{idx}); %Set the title of the axis
        %         g(gi,cond).set_color_options('chroma',3); %Set the color of the point
        %         g(gi,cond).set_color_options('map',myColors(1:3,:)/60,'n_color',3,'n_lightness',1,'lightness',60,'legend','merge'); %Set the color of the point
        g(gi,cond).set_title(sprintf('%s StepwiseRegress\n[1,2,4],coef= %.3f, p = %.3f\n[2,4],coef= %.3f, p = %.3f',groupStr{gi},mdl_124.Coefficients.Estimate(3),mdl_124.Coefficients.pValue(3),mdl_24.Coefficients.Estimate(3),mdl_24.Coefficients.pValue(3)));
        g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
        g(gi,cond).axe_property('XLim',[-0.8 0.3],'YLim',[0.5 1.1],'YTick',[0.6:0.1:1],'XTick',[-0.8:0.2:0.2]);
        g(gi,cond).set_color_options('map','brewer1');
    end
    figure('Position',[100 100 1200 550]);
    g.draw();
    saveas(gcf,fullfile(Dir.figs,['NeuVarCorr_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))
end
%% transition time point: find group average, define range = [-0.1 0.1], search individual peak time within this range
detectRange = [3.1 3.5];
clear avgPeak
for gi = 1:2
    for cond_i = 1:3
        tmp_series = mean(dat.betaCurveFrontal{gi}{cond_i});
        [pks,locs] = findpeaks(-tmp_series);
        tmpID = gndTF{1}.time(SStimeID{cond_i}(locs))<=detectRange(2) & gndTF{1}.time(SStimeID{cond_i}(locs))>detectRange(1);

        avgPeak.pks(gi,cond_i) = -pks(tmpID);
        avgPeak.locs(gi,cond_i) = gndTF{1}.time(SStimeID{cond_i}(locs(tmpID)));

    end
end

%% Delay Avg vs Probe Avg
clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontal{1},size(dat.betaAvgFrontal{1},1)*size(dat.betaAvgFrontal{1},2),1);...
    reshape(dat.betaAvgFrontalProbe{1},size(dat.betaAvgFrontalProbe{1},1)*size(dat.betaAvgFrontalProbe{1},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==1)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)]]';

stats = rm_anova2(Xa(:,1),Xa(:,4),Xa(:,2),Xa(:,3),{'Time','SS'})

mean(dat.betaAvgFrontalProbe{1})- mean(dat.betaAvgFrontal{1})+1

%% correlation between indvidual SS slope and baseline index (SS1 beha or variance)
figure('position',[100 100 700 1000]);

for gi = 1:2
    SS_slope = line_fit{gi}.params(:,1);

    [SSslopeR.r.var{gi},SSslopeR.p.var{gi}] = corr(dat.betaAvgFrontal{gi}, SS_slope,'type','Spearman');
    [SSslopeR.r.beha{gi},SSslopeR.p.beha{gi}] = corr(dat.beha{gi}, SS_slope,'type','Spearman');

    for cond_i = 1:3
        subplot(4,3,(gi-1)*3+cond_i);
        scatter(dat.betaAvgFrontal{gi}(:,cond_i),SS_slope)
        lsline
        xlabel('Delay variance')
        ylabel('Delay slope')

        title([groupStr{gi},' ',condStr{cond_i},' var'],sprintf('Rho = %.3f, p = %.3f',SSslopeR.r.var{gi}(cond_i),SSslopeR.p.var{gi}(cond_i)))
    end

    for cond_i = 1:3
        subplot(4,3,(gi-1)*3+cond_i+6);
        scatter(dat.beha{gi}(:,cond_i),SS_slope,'green')
        lsline
        xlabel('Behav')
        ylabel('Delay slope')
        title([groupStr{gi},' ',condStr{cond_i},' behav'],sprintf('Rho = %.3f, p = %.3f',SSslopeR.r.beha{gi}(cond_i),SSslopeR.p.beha{gi}(cond_i)))
    end
end
saveas(gcf,fullfile(Dir.figs,['DelayVarCorrSlope_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))

%% delay load effect topography
clear Ftest
myFigBasic;

load chanLocs_dms.mat
[is, loc] = ismember(gndTF{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.05;
figure('Name','Load-sensitive chans');
for gi = 1:2
    clear X
    for c = 1:length(gndTF{1}.label)

        load1 = squeeze(mean(mean(gndTF{gi,1}.indv(:,c,freqID.beta,timeIDdelayFtest),3),4));%subN
        load2 = squeeze(mean(mean(gndTF{gi,2}.indv(:,c,freqID.beta,timeIDdelayFtest),3),4));%subN
        load3 = squeeze(mean(mean(gndTF{gi,3}.indv(:,c,freqID.beta,timeIDdelayFtest),3),4));%subN

        X(:,1) = [load1;load2;load3];
        X(:,2) = sort(repmat([1 2 3]',sum(subs.group==gi),1));
        X(:,3) = repmat([1:sum(subs.group==gi)]',3,1);

        [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
    end
    subplot(1,2,gi);hold all
 
    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',4},'plotrad',0.53,'electrodes','off' );
    caxis([0, 6])
    h = colorbar;
    %             set(h,'ytick',[0.33:0.03:0.4])
    h.Label.String = sprintf('Ftest,p<%.3f',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];  
    title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROIdelayFtest(1),timeROIdelayFtest(2)),['N=,' num2str(sum(subs.group==gi))]);
end
           