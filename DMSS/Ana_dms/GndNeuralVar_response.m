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
IsBL2preDelay = 1;% baselined to pre response
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Young-Old'};
condStr = {'ss1','ss2','ss4'};

frontalROI = {'Fz','F1','F2','FCz','AFz'};% frontal cluster
% frontalROI = {'CPz','CP1','CP2','Pz','Cz'};%centroparietal cluster

occipROI = {'POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 12];% Hz

freq.betaFreq = [1 3];% Hz
% freq.betaFreq = [8 12];% Hz
% freq.betaFreq = [4 7];% Hz
% freq.betaFreq = [28 40];% Hz

timeROI.Post = [0.1 0.5];% in s
% timeROI.Post = [0.1 0.5];% in s
timeROI.Pre = [-0.4 -0.1];% in s
timeROI.all = [timeROI.Pre(1) timeROI.Post(2)];% in s, for curve plotting

%%
subsAll = cell(subN,3);
for sub_i = 1:subN
    subname = subs.name{sub_i};
    matName = fullfile(Dir.results,[subname,'_resp_ft_var',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

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

    cfg = [];
    cfg.operation = '(x1+x2+x3)/3';
    cfg.parameter = 'indv';
    gndTF_SSavg{gi} = ft_math(cfg, gndTF{gi,1},gndTF{gi,2},gndTF{gi,3});
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

figure('Position',[100 100 900 1200]);

for gi = 1:2
    for cond_i = 1:3
        subplot(6,3,cond_i+(gi-1)*9);
        cfg = [];
        cfg.xlim  = timeROI.all;
        cfg.ylim = [1 40];
        cfg.zlim = [-1 1];
        % cfg.colormap = bkr;
        cfg.colormap = jet;
        cfg.channel = frontalROI;
        cfg.figure = 'gca';

        ft_singleplotTFR(cfg, gndTF{gi,cond_i});

        title([groupStr{gi},condStr{cond_i}],[frontalROI{:}],'FontSize',10);
        hc = colorbar;
        hc.Label.String = 'Relative Variability(a.u.)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',timeROI.all)
        if cond_i == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time(0s=response)')
        rectangle('Position',[0 0.5 0.5 40],'LineWidth',1.2,'LineStyle','--')

        subplot(6,3,3+cond_i+(gi-1)*9);
        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        cfg.zlim = [-1 1];
        cfg.highlightchannel = frontalROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.title = sprintf('Beta(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2));
        ft_topoplotTFR(cfg, gndTF{gi,cond_i});colorbar

        subplot(6,3,6+cond_i+(gi-1)*9);
        cfg = [];
        cfg.ylim = freq.alphaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        cfg.zlim = [-1 1];
        cfg.highlightchannel = frontalROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.title = sprintf('Alpha(%d~%dHz)',freq.alphaFreq(1),freq.alphaFreq(2));

        ft_topoplotTFR(cfg, gndTF{gi,cond_i});colorbar
    end
end
saveas(gcf,fullfile(Dir.figs,['NeuVar_Resp_',[frontalROI{:}],txtCell{IsLap+1,1},num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%%  time frequency of group average
cfg = [];
cfg.xlim  = timeROI.all;
cfg.ylim = [1 40];
cfg.zlim = [-1 1];
% cfg.colormap = bkr;
cfg.colormap = jet;
cfg.channel = frontalROI;
cfg.figure = 'gca';

figure('Position',[100 100 1000 800]);
for gi = 1:3
    subplot(3,3,gi);axis square;hold all

    ft_singleplotTFR(cfg, GroupAvg{gi});

    title(groupStr{gi},'Average set sizes','FontSize',14);
    hc = colorbar;
    hc.Label.String = 'Relative Variability(a.u.)';
    hc.Label.Rotation = -90;
    hc.Label.Position = [4.2 0 0];

    xlabel('\bfTime(0s=response)')
    rectangle('Position',[0 -1 timeROI.all(2) 50],'LineWidth',1.2,'LineStyle','--')
    set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',unique([timeROI.Pre timeROI.Post]),'XLim',timeROI.all,'YLim',[0.5 40])
end

cfg = [];
cfg.ylim = freq.betaFreq;%freq
cfg.xlim = timeROI.Post;%time
cfg.zlim = [-0.4 0.4];
cfg.highlightchannel = frontalROI;
cfg.layout = 'easyCapM1';
cfg.colormap = jet;
cfg.markersymbol = '.';
cfg.highlight = 'on';
cfg.highlightsymbol = '.';
cfg.highlightsize = 18;
cfg.figure      = 'gca';

for cond_i = 1:3
    subplot(3,3,cond_i+3);axis square
    ft_topoplotTFR(cfg, gndTF_Groupdiff{cond_i});
    title(condStr{cond_i},groupStr{3});colorbar

end

for gi = 1:2
    subplot(3,3,6+gi);axis square
    ft_topoplotTFR(cfg, GroupAvg{gi});
    title(groupStr{gi});colorbar
end
saveas(gcf,fullfile(Dir.figs,['NeuVar_Resp_GroupAvg_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
saveas(gcf,fullfile(Dir.figs,['NeuVar_Resp_GroupAvg_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))

%%
freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID.all = dsearchn(gndTF{1}.time',timeROI.all');
timeID.all = timeID.all(1):timeID.all(2);

timeID.Post = dsearchn(gndTF{1}.time',timeROI.Post');
timeID.Post = timeID.Post(1):timeID.Post(2);

timeID.Pre = dsearchn(gndTF{1}.time',timeROI.Pre');
timeID.Pre = timeID.Pre(1):timeID.Pre(2);

clear chanID
[log1,chanID.frontal] = ismember(frontalROI,gndTF{1}.label);
[log2,chanID.occip] = ismember(occipROI,gndTF{1}.label);

clear dat

for gi = 1:2
    for cond_i = 1:3
        dat.betaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.betaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.betaAvgFrontalPre{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.Pre),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        dat.betaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID.all),2,"omitnan"),3,"omitnan"));
        dat.betaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.all),2,"omitnan"),3,"omitnan"));

        dat.alphaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.alphaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.alphaAvgFrontalPre{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.Pre),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        dat.alphaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID.all),2,"omitnan"),3,"omitnan"));
        dat.alphaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.all),2,"omitnan"),3,"omitnan"));
    end
end
%% mixed ANOVA: group*cond = 2*3
clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontal{1},size(dat.betaAvgFrontal{1},1)*size(dat.betaAvgFrontal{1},2),1);...
    reshape(dat.betaAvgFrontal{2},size(dat.betaAvgFrontal{2},1)*size(dat.betaAvgFrontal{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';
[SSQs.betaFrontal, DFs.betaFrontal, MSQs.betaFrontal, Fs.betaFrontal, Ps.betaFrontal]=mixed_between_within_anova(Xa);

T = subs(:,{'name','group','groupStr'});
T.load1 = [dat.betaAvgFrontal{2}(:,1);dat.betaAvgFrontal{1}(:,1)];
T.load2 = [dat.betaAvgFrontal{2}(:,2);dat.betaAvgFrontal{1}(:,2)];
T.load3 = [dat.betaAvgFrontal{2}(:,3);dat.betaAvgFrontal{1}(:,3)];
writetable(T,fullfile(Dir.ana,'TableOutput',['RespBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

Xa(:,1) = [reshape(dat.betaAvgFrontalPre{1},size(dat.betaAvgFrontalPre{1},1)*size(dat.betaAvgFrontalPre{1},2),1);...
    reshape(dat.betaAvgFrontalPre{2},size(dat.betaAvgFrontalPre{2},1)*size(dat.betaAvgFrontalPre{2},2),1)];
[SSQs.betaFrontalPre, DFs.betaFrontalPre, MSQs.betaFrontalPre, Fs.betaFrontalPre, Ps.betaFrontalPre]=mixed_between_within_anova(Xa);

Xa(:,1) = [reshape(dat.betaAvgOccip{1},size(dat.betaAvgOccip{1},1)*size(dat.betaAvgOccip{1},2),1);...
    reshape(dat.betaAvgOccip{2},size(dat.betaAvgOccip{2},1)*size(dat.betaAvgOccip{2},2),1)];
[SSQs.betaOccip, DFs.betaOccip, MSQs.betaOccip, Fs.betaOccip, Ps.betaOccip]=mixed_between_within_anova(Xa);

Xa(:,1) = [reshape(dat.alphaAvgFrontal{1},size(dat.alphaAvgFrontal{1},1)*size(dat.alphaAvgFrontal{1},2),1);...
    reshape(dat.alphaAvgFrontal{2},size(dat.alphaAvgFrontal{2},1)*size(dat.alphaAvgFrontal{2},2),1)];
[SSQs.alphaFrontal, DFs.alphaFrontal, MSQs.alphaFrontal, Fs.alphaFrontal, Ps.alphaFrontal]=mixed_between_within_anova(Xa);

Xa(:,1) = [reshape(dat.alphaAvgFrontalPre{1},size(dat.alphaAvgFrontalPre{1},1)*size(dat.alphaAvgFrontalPre{1},2),1);...
    reshape(dat.alphaAvgFrontalPre{2},size(dat.alphaAvgFrontalPre{2},1)*size(dat.alphaAvgFrontalPre{2},2),1)];
[SSQs.alphaFrontalPre, DFs.alphaFrontalPre, MSQs.alphaFrontalPre, Fs.alphaFrontalPre, Ps.alphaFrontalPre]=mixed_between_within_anova(Xa);

%% linear fit of SS
clear line_fit line_fitPre
for gi = 1:2
    for sub_i = 1:size(dat.betaAvgFrontal{gi},1)
        line_fit{gi}.params(sub_i,:) = polyfit([1 2 4]',dat.betaAvgFrontal{gi}(sub_i,:)',1);
        line_fitPre{gi}.params(sub_i,:) = polyfit([1 2 4]',dat.betaAvgFrontalPre{gi}(sub_i,:)',1);

    end
    [~,line_fit{gi}.p,~,line_fit{gi}.Tstat] = ttest(line_fit{gi}.params(:,1));
    [~,line_fitPre{gi}.p,~,line_fitPre{gi}.Tstat] = ttest(line_fitPre{gi}.params(:,1));
end
%% bar & time series -beta
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
fig =figure('Position',[100 100 800 650]);

subplot(2,3,4);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'w',HandleVisibility='off');

xcord = vertcat(hb(:).XEndPoints)';
plot(xcord(1,:),dat.betaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
for gi = 1:2
    if line_fit{gi}.p<.05
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'r','LineWidth',1.5,'HandleVisibility','off')
    else
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
    end
end
legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr,'YLim',[-inf 4])
ytickformat('%.1f')
ylabel('Variability(a.u.)')
title(sprintf('%s beta\nPost(%.1f~%.1fs)',[frontalROI{:}],timeROI.Post(1),timeROI.Post(2)))

% plot significance
if Ps.betaFrontal{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontal{gi} dat.betaAvgFrontal{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-0.7,-0.6,3);
        xposi = {[1 2],[2 3],[1 3]};
        for ss = 1:3
            if tmp_pval(ss)<=.05
                plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
                text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
            end
        end
    end
end

%-----------------inserted panel

axes('Position',[.18 .33 .1 .1],'Units','normalized')
hold all;box on;
se = std([dat.betaAvgFrontal{1};dat.betaAvgFrontal{2}])./sqrt(subN);
bar([1 2 3],mn,'stacked','BarWidth',0.5,'FaceAlpha',0.3)
errorbar([1 2 3],sum(mn),se,'k','HandleVisibility','off')
set(gca,'YLim', [0 2],'ytick',[0 2],'XTick',[1 2 3],'XTickLabel',condStr,'XLim',get(gca,'xlim')+[0.15 -0.15]);
set(gca, 'color', 'none');
if  Ps.betaFrontal{3}<.05
    plot([1 3],ones(1,2)*1.4,'k','HandleVisibility','off')
    text(2,1.5,'*','HorizontalAlignment','center','FontSize',20)
end
legend(groupStr)

%-----------------------

subplot(2,3,5);hold all;axis square
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
legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variability(a.u.)')
title(sprintf('%s beta(%.1f~%.1fs)',[occipROI{:}],timeROI.Post(1),timeROI.Post(2)))

subplot(2,3,6);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontalPre,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontalPre,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.betaAvgFrontalPre{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontalPre{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';

for gi = 1:2
    if line_fitPre{gi}.p<.05
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'r','LineWidth',1.5,'HandleVisibility','off')
    else
        errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
    end
end

legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variability(a.u.)')
title(sprintf('%s beta\nPre(%.1f~%.1fs)',[frontalROI{:}],timeROI.Pre(1),timeROI.Pre(2)))

% plot significance
if Ps.betaFrontalPre{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontalPre{gi} dat.betaAvgFrontalPre{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-0.5,-0.4,3);
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

times = gndTF{1}.time(timeID.all);
for gi = 1:2
    subplot(2,8,gi*2+4);hold all
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveFrontal{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveFrontal{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    legend(condStr,'Location','southoutside')
    ytickformat('%.1f')
    xtickangle(0)
    ylabel('Variability(a.u.)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[frontalROI{:}],freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-1.1 1.1],'YAxisLocation','right')
    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k--','HandleVisibility','off')

    subplot(2,4,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    legend(condStr,'Location','southoutside')
    ytickformat('%.1f')
    ylabel('Variability(a.u.)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-1.2 1.2])
    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
end

saveas(gca,fullfile(Dir.figs,['RespBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
fig.PaperOrientation = 'landscape';
print(fig,fullfile(Dir.figs,['RespBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r300')
%% anova F value- get fcritical from finv, permutation useing cluster test for time frequency map
clear xFs
for fi = 1:numel(gndTF{2}.freq)
    for ti = 1:length(timeID.Post)
        clear tmp_dat
        for gi = 1:2
            for cond_i = 1:3
                tmp_dat{gi,cond_i} = squeeze(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal,fi,timeID.Post(ti)),2,'omitnan')); % subs*chan*freq*time
            end
        end
        clear X
        X(:,1) = [vertcat(tmp_dat{1,:});vertcat(tmp_dat{2,:})];
        X(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
        X(:,3) = [sort(repmat([1:3]',sum(subs.group==1),1));sort(repmat([1:3]',sum(subs.group==2),1))];
        X(:,4) = [repmat([1:sum(subs.group==1)]',3,1);repmat([1:sum(subs.group==2)]',3,1)+sum(subs.group==1)];

        [~, ~, ~, xFs{fi,ti}, ~]=mixed_between_within_anova(X);
    end
end

clear xFsMat
for fi = 1:numel(gndTF{2}.freq)
    for ti = 1:length(timeID.Post)
        xFsMat.inter(fi,ti) = xFs{fi,ti}{4};
        xFsMat.post(fi,ti) = xFs{fi,ti}{3};
    end
end

threshP_Finter = 0.05;

Fcritical = finv(1-threshP_Finter,2,78);% for interaction

figure('position',[200 200 800 500]);

for i = 1:2
    subplot(1,2,i);hold all;axis square
    if i ==1
        tmpdat = xFsMat.inter>Fcritical;
        title('Interaction effect binarized Fmap')
    else
        tmpdat = xFsMat.post>Fcritical;
        title('Load main effect binarized Fmap')
    end
    imagesc(gndTF{1}.time(timeID.Post),gndTF{1}.freq,tmpdat)
    % im.AlphaData = xFsMat.post>Fcritical;

    xlabel('Time(0s=response)');ylabel('Frequency(Hz)')
    set(gca,'ylim',[1 40],'xlim',[0 0.5])

    % caxis([Fcritical, 6])
    h = colorbar;
    h.Label.String = sprintf('Ftest(p<%.3f)',threshP_Finter);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
end
saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTimeFreqFinteraction_',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

%% bar & time series -alpha
fig =figure('Position',[100 100 800 650]);

subplot(2,3,4);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
plot(xcord(1,:),dat.alphaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

for gi = 1:2
    errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
end
legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variability(a.u.)')
title(sprintf('%s alpha\nPost(%.1f~%.1fs)',[frontalROI{:}],timeROI.Post(1),timeROI.Post(2)))
% % plot significance
% if Ps.alphaFrontal{3}<.05
%     for gi = 1:2
%         tmpdata = [dat.alphaAvgFrontal{gi} dat.alphaAvgFrontal{gi}(:,1)];
%         tmpdata = diff(tmpdata,1,2);
%         [~,tmp_pval]= ttest(tmpdata);
%         sigH = linspace(-0.8,-0.6,3);
%         xposi = {[1 2],[2 3],[1 3]};
%         for ss = 1:3
%             if tmp_pval(ss)<=.05
%                 plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
%                 text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
%             end
%         end
%     end
% end

subplot(2,3,5);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgOccip,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgOccip,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.alphaAvgOccip{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgOccip{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variability(a.u.)')
title(sprintf('%s alpha(%.1f~%.1fs)',[occipROI{:}],timeROI.Post(1),timeROI.Post(2)))

subplot(2,3,6);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgFrontalPre,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgFrontalPre,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.alphaAvgFrontalPre{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgFrontalPre{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';

for gi = 1:2
    errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
end

legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variability(a.u.)')
title(sprintf('%s alpha\nPre(%.1f~%.1fs)',[frontalROI{:}],timeROI.Pre(1),timeROI.Pre(2)))

% plot significance
if Ps.alphaFrontalPre{3}<.05
    for gi = 1:2
        tmpdata = [dat.alphaAvgFrontalPre{gi} dat.alphaAvgFrontalPre{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-0.5,-0.4,3);
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
times = gndTF{1}.time(timeID.all);
for gi = 1:2
    subplot(2,8,gi*2+4);hold all
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveFrontal{gi}(:,cond_i,:)));
        se = squeeze(std(dat.alphaCurveFrontal{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    legend(condStr,'Location','southoutside')
    ytickformat('%.1f')
    xtickangle(0)
    ylabel('Variability(a.u.)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[frontalROI{:}],freq.alphaFreq(1),freq.alphaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-1.2 1.2],'YAxisLocation','right')
    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k--','HandleVisibility','off')

    subplot(2,4,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.alphaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    legend(condStr,'Location','southoutside')
    ytickformat('%.1f')
    ylabel('Variability(a.u.)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.alphaFreq(1),freq.alphaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-1.2 1.2])
    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
end

saveas(gca,fullfile(Dir.figs,['RespAlphaVariance_',num2str(freq.alphaFreq(1)),'~',num2str(freq.alphaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);
behaStr = {'RT','acc'};

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',2);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

for idx = 2%1:2
    if idx == 2
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));
    else
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_rt','s2_rt','s3_rt'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_rt','s2_rt','s3_rt'}));
    end

    clear g

    Y = [mean(dat.beha{1},2);mean(dat.beha{2},2)];%old=2,young= 1
    X = [mean(dat.betaAvgFrontal{1},2);mean(dat.betaAvgFrontal{2},2)];
    condArray = flipud(subs.group);

    %  a = X(condArray==2);
    %     b = Y(condArray==2);
    %     [brob,stats] = robustfit(a,b);

    [r(1),pval(1)] = corr(X(condArray==1),Y(condArray==1),'type','Spearman');%young
    [r(2),pval(2)] = corr(X(condArray==2),Y(condArray==2),'type','Spearman');%old

    %     [r(1),pval(1)] = corr(X(condArray==1),Y(condArray==1),'type','Pearson');%young
    %     [r(2),pval(2)] = corr(X(condArray==2),Y(condArray==2),'type','Pearson');%old

    g(1) = gramm('x',X,'y',Y,'color',groupStr(condArray)); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
    g(1).geom_point(); %Draw a scatter plot
    g(1).stat_glm(); %Draw lines and confidence intervals based on scatter plots
    g(1).set_names('x',sprintf('Beta variance(%d~%dHz)\n%s',freq.betaFreq(1),freq.betaFreq(2),[frontalROI{:}]),'y',behaStr{idx}); %Set the title of the axis
    g(1).set_color_options('map',myColors(1:2,:),'n_lightness',1); %Set the color of the point
    g(1).set_title(sprintf('collapse set sizes\nOld, Rho = %.3f, p = %.3f\nYoung, Rho = %.3f, p = %.3f',r(2),pval(2),r(1),pval(1)));
    g(1).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    g(1).axe_property('XLim',[-0.5 2.5]);

    X = [mean(dat.alphaAvgFrontal{1},2);mean(dat.alphaAvgFrontal{2},2)];

    [r(1),pval(1)] = corr(X(condArray==1),Y(condArray==1),'type','Spearman');%young
    [r(2),pval(2)] = corr(X(condArray==2),Y(condArray==2),'type','Spearman');%old
    g(2) = gramm('x',X,'y',Y,'color',groupStr(condArray)); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
    g(2).geom_point(); %Draw a scatter plot
    g(2).stat_glm(); %Draw lines and confidence intervals based on scatter plots
    g(2).set_names('x',sprintf('Alpha variance(%d~%dHz)\n%s',freq.alphaFreq(1),freq.alphaFreq(2),[frontalROI{:}]),'y',behaStr{idx}); %Set the title of the axis
    g(2).set_color_options('map',myColors(1:2,:),'n_lightness',1); %Set the color of the point
    g(2).set_title(sprintf('collapse set sizes\nOld, Rho = %.3f, p = %.3f\nYoung, Rho = %.3f, p = %.3f',r(2),pval(2),r(1),pval(1)));
    g(2).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    g(2).axe_property('XLim',[-1 3]);
    figure('Position',[100 100 500 260]);
    g.draw();
    saveas(gcf,fullfile(Dir.figs,['RespNeuVarCorr_',behaStr{idx},num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
end

myFigBasic;

load chanLocs_dms.mat
[is, loc] = ismember(gndTF{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));

myColors = jet;

% myColors = crameri('vik',256);
% myColors = myColors(256/2:256,:);

%% glme regression model: acc = age-group*beta variability
tmp = vertcat(dat.betaAvgFrontal{:});
betaArray = [tmp(:,1);tmp(:,2);tmp(:,3)];
groupArray = flipud([subs.group;subs.group;subs.group]);
subArray = repmat([1:subN],1,3)';
loadArray = [ones(subN,1)*1;ones(subN,1)*2;ones(subN,1)*3];
tmp = vertcat(dat.beha{:});
accArray = [tmp(:,1);tmp(:,2);tmp(:,3)];

tbl = table(accArray,betaArray, subArray, loadArray,groupArray,'VariableNames',{'beha','var','subj','ss','group'});
tbl.subj = categorical(tbl.subj);
tbl.group = categorical(tbl.group);
tbl.ss = categorical(tbl.ss);
X = [betaArray,loadArray,groupArray];
y = accArray;

glme = fitglme(tbl,'beha ~ var*group+(1|ss)');
% glme = fitglme(tbl,'beha ~ var*group+(1|ss)','Distribution','Gamma');
t = dataset2table(glme.anova);
t.fomula{1} = char(glme.Formula); 
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova'])

t = dataset2table(glme.Coefficients);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients'])

tbl(double(tbl.subj)==35,:)=[]; % remove the outlier for frontal beta
glme = fitglme(tbl,'beha ~ var*group+(1|ss)');
t = dataset2table(glme.anova);

full_model = fitglme(tbl,'beha ~ var*group+(1|ss)');
reduced_model = fitglme(tbl,'beha ~ var+group+(1|ss)');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood);

% t.fomula{1} = char(glme.Formula); 
% writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaVarglme_wo_outlier',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%     'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova'])
% 
% t = dataset2table(glme.Coefficients);
% writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaVarglme_wo_outlier',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%     'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients'])

%% normalize accuracy distribution: accuracy was a binary variable, averaging across trials make a continous variable ranging [0 1]
% our sample's avg acc has a beta distribution; therefore, we fit a beta
% distribution to estimate a & b, then use the predicted value to run glme

tbl.beha_Bnorm = tbl.beha;
params = betafit(tbl.beha);
tbl.beha_Bnorm = betainv(tbl.beha, params(1), params(2));
writetable(tbl,fullfile(Dir.ana,'TableOutput',['RespBetaVariance_Regress_YoungOld_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

% figure;histogram(tbl.beha);hold on;histogram(tbl.beha_Bnorm);
% xlabel('Accuracy');ylabel('Density');
% saveas(gca,fullfile(Dir.figs,'BetaNormalizationAvgACC.png'))

glme = fitglme(tbl,'beha_Bnorm ~ var*group +(1|ss)');

t = dataset2table(glme.anova);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova_Bnorm'])

t = dataset2table(glme.Coefficients);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coef_Bnorm'])

full_model = fitglme(tbl,'beha_Bnorm ~ var*group +(1|ss)');
reduced_model = fitglme(tbl,'beha_Bnorm ~ var+group +(1|ss)');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood)

%%
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',5);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
myColors = [237 28 50;0 166 69;25 119 200 ]./255;

for idx = 2%1:2
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
        %             [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Pearson');

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
        mdl_124 =  fitglm(X,Y);
        mdl_124_0 = fitglm(X(:,1),Y);

        T = subs(subs.group ==gi,{'name','group','groupStr'});
        T = [T;T;T];% repeat three times.
        T.respVar = X(:,2);T.ss = X(:,1);T.beha = Y;
        writetable(T,fullfile(Dir.ana,'TableOutput',['RespBetaVariance_Regress_',groupStr{gi},behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))
%         tbl{gi} = T;
        T.group = categorical(T.group);

        clear X Y
        Y = reshape(dat.beha{gi}(:,[2 3]),size(dat.beha{gi},1)*2,1);
        X(:,1) = sort(repmat([2 4]',size(dat.betaAvgFrontal{gi},1),1));
        X(:,2) = reshape(dat.betaAvgFrontal{gi}(:,[2 3]),size(dat.betaAvgFrontal{gi},1)*2,1);

        mdl_24 =  fitglm(X,Y);
        mdl_24_0 = fitglm(X(:,1),Y);
        
        if gi == 2
            clear X Y
            Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
            X(:,1) = sort(repmat([1 2 4]',size(dat.betaAvgFrontal{gi},1),1));
            X(:,2) = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);
            X([15 36 57],:) = [];
            Y([15 36 57],:) = [];
            mdl_124_wo =  fitglm(X,Y);
            mdl_124_wo_0 =  fitglm(X(:,1),Y);

            clear X Y
            Y = reshape(dat.beha{gi}(:,[2 3]),size(dat.beha{gi},1)*2,1);
            X(:,1) = sort(repmat([2 4]',size(dat.betaAvgFrontal{gi},1),1));
            X(:,2) = reshape(dat.betaAvgFrontal{gi}(:,[2 3]),size(dat.betaAvgFrontal{gi},1)*2,1);
            X([15 36],:) = [];
            Y([15 36],:) = [];
            mdl_24_wo =  fitglm(X,Y);
            mdl_24_wo_0 =  fitglm(X(:,1),Y);
        end
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
        g(gi,cond).set_color_options('map',myColors(1:3,:)/60,'n_color',3,'n_lightness',1,'lightness',60,'legend','merge'); %Set the color of the point

%         g(gi,cond).set_title(sprintf('%s\n[1,2,4]wo,B= %.3f, p = %.3f, Rsq= %.3f\n[2,4]wo,B = %.3f, p = %.3f, Rsq= %.3f',groupStr{gi},mdl_124_wo.Coefficients.Estimate(3),mdl_124_wo.Coefficients.pValue(3),mdl_124_wo.Rsquared.Ordinary-mdl_124_wo_0.Rsquared.Ordinary,mdl_24_wo.Coefficients.Estimate(3),mdl_24_wo.Coefficients.pValue(3),mdl_24_wo.Rsquared.Ordinary-mdl_24_wo_0.Rsquared.Ordinary));

        g(gi,cond).set_title(sprintf('%s\n[1,2,4],B= %.3f, p = %.3f, Rsq= %.3f\n[2,4],B = %.3f, p = %.3f, Rsq= %.3f',groupStr{gi},mdl_124.Coefficients.Estimate(3),mdl_124.Coefficients.pValue(3),mdl_124.Rsquared.Ordinary-mdl_124_0.Rsquared.Ordinary,mdl_24.Coefficients.Estimate(3),mdl_24.Coefficients.pValue(3),mdl_24.Rsquared.Ordinary-mdl_24_0.Rsquared.Ordinary));

        g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times

        if idx == 2
            g(gi,cond).axe_property('YLim',[0.6 1],'YTick',[0.6:0.1:1]);
        else
            g(gi,cond).axe_property('YLim',[0.4 1.4],'YTick',[0.4:0.2:2]);
        end
        if gi ==1
            g(gi,cond).axe_property('XLim',[-1 2],'XTick',[-1:1:4]);
        else
            g(gi,cond).axe_property('XLim',[-1 3],'XTick',[-1:1:4]);
        end
        g(gi,cond).set_color_options('map',myColors(1:3,:),'n_lightness',1);
        g(gi,cond).set_point_options('base_size',2);
        g(gi,cond).set_line_options('base_size',0.5);
    end
    fig = figure('Position',[100 100 1000 500]);
    g.draw();
    saveas(gcf,fullfile(Dir.figs,['RespVarGLM_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
    fig.PaperOrientation = 'landscape';
    print(fig,fullfile(Dir.figs,['RespVarGLM_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']),'-dpdf','-r300','-bestfit')
end
return
%% beta rebound relative to pre-response, average cross set sizes
figure('Name','Ftest chans','Position',[100 200 1300 300]);
threshP = 0.001;
mkSize = 3;
groupLimit = {[0,80],[0 80]};
myColors = flipud(hot);
clear Ftest

for gi = 1:2
    for c = 1:length(gndTF{1}.label)

        pre = squeeze(mean(mean(gndTF_SSavg{gi}.indv(:,c,freqID.beta,timeID.Pre),3),4));%subN
        post = squeeze(mean(mean(gndTF_SSavg{gi}.indv(:,c,freqID.beta,timeID.Post),3),4));%subN
        clear X

        X(:,1) = [pre;post];
        X(:,[3 2]) = fullfact([sum(subs.group==gi) 2]);

        [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
    end
    subplot(1,5,gi);hold all

    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors,'style','contour');
    caxis(groupLimit{gi})
    h = colorbar;
    h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.Post(1),timeROI.Post(2)),['N=' num2str(sum(subs.group==gi))]);
end
orient landscape
print(gcf,fullfile(Dir.figs,['RespVarTopoMixedANOVA22_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']),'-dtiff')

%% 2*2 mixed anova 2(Pre vs Post)*2(Young vs Old)
threshP = 0.05;
clear Ftest
for c = 1:length(gndTF{1}.label) % 2(Pre vs Post)*2(Young vs Old)
    pre = squeeze(mean(mean(gndTF_SSavg{1}.indv(:,c,freqID.beta,timeID.Pre),3),4));%Young
    tmp = squeeze(mean(mean(gndTF_SSavg{2}.indv(:,c,freqID.beta,timeID.Pre),3),4));%Old
    pre = [pre;tmp];

    post = squeeze(mean(mean(gndTF_SSavg{1}.indv(:,c,freqID.beta,timeID.Post),3),4));%young
    tmp = squeeze(mean(mean(gndTF_SSavg{2}.indv(:,c,freqID.beta,timeID.Post),3),4));%old
    post = [post;tmp];

    clear X
    X(:,1) = [pre;post];
    X(:,2) = flipud([subs.group;subs.group]);
    X(:,[4 3]) = fullfact([subN 2]);

    [~, ~, ~, Ftest.F{c},Ftest.p{c}]=mixed_between_within_anova(X);
end

clear xFsMat
for c = 1:length(gndTF{1}.label)
    xFsMat.inter(c) = Ftest.F{c}{4};
    xFsMat.post(c) = Ftest.F{c}{3};
    xFsMat.group(c) = Ftest.F{c}{1};

    xFsMat.inter_p(c) = Ftest.p{c}{4};
    xFsMat.post_p(c) = Ftest.p{c}{3};
    xFsMat.group_p(c) = Ftest.p{c}{1};
end

subplot(1,5,5);hold all % 2*2 interaction
topoplot(xFsMat.inter,chanloc,'emarker2',{find(xFsMat.inter_p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
caxis([0, 6])
h = colorbar;
h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
h.Label.Rotation = -90;
h.Label.Position = h.Label.Position + [1 0 0];
title(sprintf('interaction\n%.1f~%.1fs', timeROI.Post(1),timeROI.Post(2)));

subplot(1,5,3);hold all
topoplot(xFsMat.post,chanloc,'emarker2',{find(xFsMat.post_p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
caxis([0, 80])
h = colorbar;
h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
h.Label.Rotation = -90;
h.Label.Position = h.Label.Position + [1 0 0];
title(sprintf('Pre vs Post main\n%.1f~%.1fs', timeROI.Post(1),timeROI.Post(2)));

subplot(1,5,4);hold all
topoplot(xFsMat.group,chanloc,'emarker2',{find(xFsMat.group_p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
caxis([0, 6])
h = colorbar;
h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
h.Label.Rotation = -90;
h.Label.Position = h.Label.Position + [1 0 0];
title(sprintf('Young vs Old\n%.1f~%.1fs', timeROI.Post(1),timeROI.Post(2)));

saveas(gcf,fullfile(Dir.figs,['RespVarTopoMixedANOVA22_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
%% spatial cluster permutation on pre vs post, separate age group
threshP = 0.001;
for gi = 1:2
    ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));

    cfg= [];
    cfg.latency = 0;
    cfg.keepindividual = 'yes';
    tmpsubs = subs.group==gi;
    cfg.channel = tfDat{1}.label;
    ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{tmpsubs,1});

    cfg = [];
    cfg.layout = 'easycapM11.mat';
    layout = ft_prepare_layout(cfg);

    cfg = [];
    cfg.layout = layout;
    cfg.method = 'triangulation';
    neighbours = ft_prepare_neighbours(cfg);

    ft_post = ft_temp;
    ft_post.individual = squeeze(mean(mean(gndTF_SSavg{gi}.indv(:,:,freqID.beta,timeID.Post),3),4));%subN*chan

    ft_pre = ft_temp;
    ft_pre.individual = squeeze(mean(mean(gndTF_SSavg{gi}.indv(:,:,freqID.beta,timeID.Pre),3),4));%subN*chan

    cfg                  = [];
    cfg.parameter        = 'individual';
    cfg.neighbours       = neighbours;
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusterstatistic = 'maxsum';
    cfg.clusteralpha     = threshP;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = threshP;
    cfg.minnbchan        = 2;%% minimal number of neighbouring channels

    cfg.numrandomization = 10000;

    Nsubj  = sum(subs.group==gi);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2; % number or list with indices indicating the independent variable(s)

    chanStat{gi}              = ft_timelockstatistics(cfg, ft_post, ft_pre);
end

%%  make a plot

for gi = 1:2
    figure;
    cfg = [];
    cfg.layout = 'easycapM11.mat';
    cfg.parameter='stat';
    cfg.contournum = 0;
    cfg.zlim = [0 10];
    cfg.comment = 'no';
    if sum(chanStat{gi}.mask)>0

        cfg.highlightsymbolseries = ['*','*','.','.','.'];
        cfg.markersymbol = '.';
        cfg.alpha = threshP;
        cfg.subplotsize = [1 1];
        ft_clusterplot(cfg, chanStat{gi} );
    else
        cfg.marker             = 'off';
        ft_topoplotER(cfg, chanStat{gi})
    end

    colormap(jet);
    h = colorbar;
    h.Ticks = [0:5:10];  
    h.Label.String = sprintf('T-value(cluster permutation with p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(groupStr{gi})
    saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoPrePostpermutation_',groupStr{gi},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
    saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoPrePostpermutation_',groupStr{gi},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))
end

%% post-response beta rebound, 2*3 mixed anova, group*set sizes
figure('Name','Ftest chans','Position',[100 200 800 300]);
threshP = 0.05;
mkSize = 3;
groupLimit = {[0,80],[0 80]};

clear Ftest

for c = 1:length(gndTF{1}.label)
    tmp = [];
    for i = 1:3% set sizes
        g1 = squeeze(mean(mean(gndTF{1,i}.indv(:,c,freqID.beta,timeID.Post),3),4));%subN
        g2 = squeeze(mean(mean(gndTF{2,i}.indv(:,c,freqID.beta,timeID.Post),3),4));%subN

        tmp = [tmp;g1;g2];
    end
    clear X
    X(:,1)  =tmp;
    X(:,[4 3]) = fullfact([subN 3]);
    X(:,2) = flipud(repmat(subs.group,3,1));

    [~, ~, ~, Ftest.F{c}, Ftest.p{c}]=mixed_between_within_anova(X);

    Ftest.interF(c) = Ftest.F{c}{4};
    Ftest.loadF(c) = Ftest.F{c}{3};
    Ftest.groupF(c) = Ftest.F{c}{1};

    Ftest.inter_p(c) = Ftest.p{c}{4};
    Ftest.load_p(c) = Ftest.p{c}{3};
    Ftest.group_p(c) = Ftest.p{c}{1};
end

for i = 1:3
    subplot(1,3,i);hold all
    switch i
        case 1
            topoplot(Ftest.interF,chanloc,'emarker2',{find(Ftest.inter_p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off');
            title(sprintf('Interaction\n%.1f~%.1fs', timeROI.Post(1),timeROI.Post(2)));
        case 2
            topoplot(Ftest.loadF,chanloc,'emarker2',{find(Ftest.load_p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off');
            title(sprintf('Load main\n%.1f~%.1fs', timeROI.Post(1),timeROI.Post(2)));
        case 3
            topoplot(Ftest.groupF,chanloc,'emarker2',{find(Ftest.group_p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off');
            title(sprintf('Group main\n%.1f~%.1fs', timeROI.Post(1),timeROI.Post(2)));
    end

    caxis([0 6])
    h = colorbar;
    h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
end
saveas(gcf,fullfile(Dir.figs,['RespVarTopoMixedANOVA23_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%% median split young and old based on their acc

beha.accAvg = mean([beha.s1_acc,beha.s2_acc, beha.s3_acc],2);
figure('position',[200 200 600 600]);
splitStr = {'High-performer','Low-performer'};
threshP = .01;
myColormap = jet;
clear chanTtest
for gi = 1:2
    tmpsubs = find(subs.group == gi);
    tmpACC = beha.accAvg(tmpsubs);

    for b = 1:2
        if b == 1
            tmp_upper = tmpACC>median(tmpACC);
        else
            tmp_upper = tmpACC<=median(tmpACC);
        end

        subplot(2,2,(gi-1)*2+b)

        for c = 1:length(gndTF{1}.label)
            pre = squeeze(mean(mean(gndTF{gi,1}.indv(tmp_upper,c,freqID.beta,timeID.Pre),3),4));%subN
            post = squeeze(mean(mean(gndTF{gi,2}.indv(tmp_upper,c,freqID.beta,timeID.Post),3),4));%subN
            [~,chanTtest.p(gi,c),~,s] = ttest(post,pre);
            chanTtest.T(gi,c) = s.tstat;
            chanTtest.data(gi,c) = mean(post);
        end

        topoplot(chanTtest.data(gi,:),chanloc,'emarker2',{find(chanTtest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColormap,'style','map');
        caxis([0 1])
        h = colorbar;
        h.Label.String = sprintf('Beta variability(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];


        %         topoplot(chanTtest.T(gi,:),chanloc,'emarker2',{find(chanTtest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off');
        %         caxis([0 7])
        %         h = colorbar;
        %         h.Label.String = sprintf('Ttest(p<%.3f)',threshP);
        %         h.Label.Rotation = -90;
        %         h.Label.Position = h.Label.Position + [1 0 0];
        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.Post(1),timeROI.Post(2)),[splitStr{b},' N=' num2str(sum(tmp_upper))]);
    end
end
saveas(gcf,fullfile(Dir.figs,['RespVarTopoTtestSplitACC_',num2str(threshP),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
%%
threshP = 0.001;
for gi = 1:2
    tmpsubs = find(subs.group == gi);
    tmpACC = beha.accAvg(tmpsubs);

    for b = 1:2
        if b == 1
            tmp_upper = tmpACC>median(tmpACC);
        else
            tmp_upper = tmpACC<=median(tmpACC);
        end

        ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));
        cfg= [];
        cfg.latency = 0;
        cfg.keepindividual = 'yes';
        cfg.channel = tfDat{1}.label;
        ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{tmp_upper,1});

        cfg = [];
        cfg.layout = 'easycapM11.mat';
        layout = ft_prepare_layout(cfg);

        cfg = [];
        cfg.layout = layout;
        cfg.method = 'triangulation';
        neighbours = ft_prepare_neighbours(cfg);

        ft_post = ft_temp;
        ft_post.individual = squeeze(mean(mean(gndTF{gi,1}.indv(tmp_upper,:,freqID.beta,timeID.Post),3),4));%subN

        ft_pre = ft_temp;
        ft_pre.individual = squeeze(mean(mean(gndTF{gi,2}.indv(tmp_upper,:,freqID.beta,timeID.Pre),3),4));%subN

        cfg                  = [];
        %         cfg.channel = 'All';
        cfg.parameter = 'individual';
        %         cfg_neighb.method = 'distance';
        cfg.neighbours       = neighbours;
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'ft_statfun_depsamplesT';
        cfg.correctm         = 'cluster';

        cfg.clusteralpha     = threshP;
        cfg.clusterstatistic = 'maxsum';
        cfg.minnbchan        = 2;%% minimal number of neighbouring channels
        cfg.tail             = 0;
        cfg.clustertail      = 0;
        cfg.alpha            = threshP;
        cfg.numrandomization = 10000;

        nsubj                     = sum(tmp_upper);
        design                    = zeros(2,2*nsubj);
        design(1,1:nsubj)         = 1;
        design(1,nsubj+1:2*nsubj) = 2;
        design(2,1:nsubj)         = 1:nsubj;
        design(2,nsubj+1:2*nsubj) = 1:nsubj;

        cfg.design   = design; % design matrix
        cfg.ivar     = 1;      % the 1st row codes the independent variable (sedation level)
        cfg.uvar     = 2;      % the 2nd row codes the unit of observation (subject)

        chanStat{gi,b}              = ft_timelockstatistics(cfg, ft_post, ft_pre);
    end
end
%%  make a plot

for gi = 1:2
    for b = 1:2
        figure('position',[200 200 300 300]);
        cfg = [];
        cfg.layout = 'easycapM11.mat';
        cfg.parameter='stat';
        cfg.contournum = 0;
        cfg.zlim = [0 10];
        cfg.comment = 'no';
        if sum(chanStat{gi,b}.mask)>0

            cfg.highlightsymbolseries = ['*','*','.','.','.'];
            cfg.markersymbol = '.';
            cfg.alpha = threshP;
            cfg.subplotsize = [1 1];
            ft_clusterplot(cfg, chanStat{gi,b} );
        else
            %                 cfg.marker             = 'off';
            cfg.markersize=0.000001;
            ft_topoplotER(cfg, chanStat{gi,b})
        end

        colormap(hot);
        h = colorbar;
        h.Label.String = sprintf('T-value(cluster permutation with p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];

        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.Post(1),timeROI.Post(2)),splitStr{b});
        saveas(gcf,fullfile(Dir.figs,['RespVarTopoTtestSplitACC_clusterPerm',groupStr{gi},num2str(b),num2str(threshP),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
        print(gcf,fullfile(Dir.figs,['RespVarTopoTtestSplitACC_clusterPerm',groupStr{gi},num2str(b),num2str(threshP),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']),'-dpdf','-r300')
    end
end
