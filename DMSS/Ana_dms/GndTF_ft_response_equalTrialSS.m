clear;
rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 1;
IsBL2preDelay = 1;% baselined to preResp
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
groupStr = {'Young','Old','Young-Old'};
condStr = {'ss1','ss2','ss4'};

% frontalROI = {'CPz','CP1','CP2','Pz','Cz'};
frontalROI = {'Fz','F1','F2','FCz','AFz'};

occipROI = {'POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz

timeROI.all = [-0.4 0.5];% in s, for curve plotting
% timeROI.all = [-0.4 0.95];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s
% timeROI.Post = [0 0.95];% in s
timeROI.Pre = [-0.4 0];% in s
%%
subsAll = cell(subN,3);

for sub_i = 1:subN
    subname = subs.name{sub_i};

    matName = fullfile(Dir.results,[subname,'_resp_ft_equalizedTrials',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
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
    gndTFdiff{gi,1} = ft_math(cfg, gndTF{gi,2},gndTF{gi,1}); %set size difference
    gndTFdiff{gi,2} = ft_math(cfg, gndTF{gi,3},gndTF{gi,1});
    gndTFdiff{gi,3} = ft_math(cfg, gndTF{gi,3},gndTF{gi,2});
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
figure('Position',[100 100 900 800]);

for gi = 1:2 
    for cond_i = 1:3
        subplot(4,3,cond_i+(gi-1)*6);

        cfg = [];
        cfg.xlim  = timeROI.all;
        cfg.ylim = [1 40];
        cfg.zlim = [-3 3];
        % cfg.colormap = bkr;
        cfg.colormap = jet;

        cfg.channel = frontalROI;
        cfg.figure = 'gca';

        ft_singleplotTFR(cfg, gndTF{gi,cond_i});

        title(condStr{cond_i},groupStr{gi},'FontSize',14);
        hc = colorbar;
        hc.Label.String = 'Power(dB)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',[0 3])
        if cond_i == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time(0s=response)')
        rectangle('Position',[0 0.5 0.5 40],'LineWidth',1.2,'LineStyle','--')

        subplot(4,3,3+cond_i+(gi-1)*6);

        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        cfg.zlim = [-3 3];
        cfg.highlightchannel = frontalROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        % cfg.marker = 'labels';
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';

        ft_topoplotTFR(cfg, gndTF{gi,cond_i});
        colorbar     

    end
    saveas(gcf,fullfile(Dir.figs,['Resp_power_equal_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4}, '.png']))
end

%%  time frequency of group average
cfg = [];
cfg.xlim  = timeROI.all;
cfg.ylim = [1 40];
cfg.zlim = [-1.5 1.5];
% cfg.colormap = bkr;
cfg.colormap = jet;
cfg.channel = frontalROI;
cfg.figure = 'gca';

figure('Position',[100 100 1000 500]);
for gi = 1:3
    subplot(2,3,gi);    axis square;hold all

    ft_singleplotTFR(cfg, GroupAvg{gi});

    title(groupStr{gi},'Average set sizes','FontSize',14);
    hc = colorbar;
    hc.Label.String = 'Power(dB)';
    hc.Label.Rotation = -90;
    hc.Label.Position = [4.2 0 0];

    xlabel('Time(0s=response)')
    rectangle('Position',[0 -0.5 timeROI.all(2) 50],'LineWidth',1.2,'LineStyle','--')
    set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',timeROI.all,'YLim',[0.5 40])
end

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

for cond_i = 1:3
    subplot(2,3,cond_i+3);axis square
    ft_topoplotTFR(cfg, gndTF_Groupdiff{cond_i});
    title(condStr{cond_i},groupStr{3});colorbar
end

saveas(gcf,fullfile(Dir.figs,['Resp_power_equal_GroupAvg_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4}, '.png']))

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

clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontalPre{1},size(dat.betaAvgFrontalPre{1},1)*size(dat.betaAvgFrontalPre{1},2),1);...
    reshape(dat.betaAvgFrontalPre{2},size(dat.betaAvgFrontalPre{2},1)*size(dat.betaAvgFrontalPre{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.betaFrontalPre, DFs.betaFrontalPre, MSQs.betaFrontalPre, Fs.betaFrontalPre, Ps.betaFrontalPre]=mixed_between_within_anova(Xa);

clear Xa
Xa(:,1) = [reshape(dat.betaAvgOccip{1},size(dat.betaAvgOccip{1},1)*size(dat.betaAvgOccip{1},2),1);...
    reshape(dat.betaAvgOccip{2},size(dat.betaAvgOccip{2},1)*size(dat.betaAvgOccip{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

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
legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Power(dB)')
title(sprintf('%s beta\nPost(%.1f~%.1fs)',[frontalROI{:}],timeROI.Post(1),timeROI.Post(2)))

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
ylabel('Power(dB)')
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
ylabel('Power(dB)')
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
    ylabel('Power(dB)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[frontalROI{:}],freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-0.5 3],'YAxisLocation','right')
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
    ylabel('Power(dB)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-1 3])
    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
end
saveas(gca,fullfile(Dir.figs,['RespBetaPower_equal_',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
%%
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);
behaStr = {'RT','acc'};

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',2);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

for idx = 1:2
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

    [r(1),pval(1)] = corr(X(condArray==1),Y(condArray==1),'type','Spearman');%young
    [r(2),pval(2)] = corr(X(condArray==2),Y(condArray==2),'type','Spearman');%old
    g(1) = gramm('x',X,'y',Y,'color',groupStr(condArray)); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
    g(1).geom_point(); %Draw a scatter plot
    g(1).stat_glm(); %Draw lines and confidence intervals based on scatter plots
    g(1).set_names('x',sprintf('Beta power(%d~%dHz)\n%s',freq.betaFreq(1),freq.betaFreq(2),[frontalROI{:}]),'y',behaStr{idx}); %Set the title of the axis
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
    g(2).set_names('x',sprintf('Alpha power(%d~%dHz)\n%s',freq.alphaFreq(1),freq.alphaFreq(2),[frontalROI{:}]),'y',behaStr{idx}); %Set the title of the axis
    g(2).set_color_options('map',myColors(1:2,:),'n_lightness',1); %Set the color of the point
    g(2).set_title(sprintf('collapse set sizes\nOld, Rho = %.3f, p = %.3f\nYoung, Rho = %.3f, p = %.3f',r(2),pval(2),r(1),pval(1)));
    g(2).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    g(2).axe_property('XLim',[-1 3]);
    figure('Position',[100 100 500 260]);
    g.draw();
      saveas(gcf,fullfile(Dir.figs,['RespBetaPowerBehavCorr_equal_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
    saveas(gcf,fullfile(Dir.figs,['RespBetaPowerBehavCorr_equal_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))
end

%% resp topography

myFigBasic;

load chanLocs_dms.mat
[is, loc] = ismember(gndTF{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.05;
mkSize = 1;
myColors = jet;

% myColors = crameri('vik',256);
% myColors = myColors(256/2:256,:);

figure('Name','Load-sensitive chans','Position',[100 200 1000 600]);
timeROI.bins = {[-0.4 0];[0 0.5]};

for t = 1:length(timeROI.bins) % loop time roi
    tmpID = dsearchn(gndTF{1}.time',timeROI.bins{t}');
    tmpID = tmpID(1):tmpID(2);

    for gi = 1:2
        clear X
        for c = 1:length(gndTF{1}.label)

            load1 = squeeze(mean(mean(gndTF{gi,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
            load2 = squeeze(mean(mean(gndTF{gi,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
            load3 = squeeze(mean(mean(gndTF{gi,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN

            X(:,1) = [load1;load2;load3];
            X(:,2) = sort(repmat([1 2 3]',sum(subs.group==gi),1));
            X(:,3) = repmat([1:sum(subs.group==gi)]',3,1);

            [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
        end
        subplot(3,length(timeROI.bins),(gi-1)*length(timeROI.bins)+t);hold all

        topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off' );
        caxis([0, 10])
        h = colorbar;
        h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.bins{t}(1),timeROI.bins{t}(2)),['N=' num2str(sum(subs.group==gi))]);
    end

    clear X
    for c = 1:length(gndTF{1}.label)

        load1 = squeeze(mean(mean(gndTF{1,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load1 = [load1;tmp];

        load2 = squeeze(mean(mean(gndTF{1,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load2 = [load2;tmp];

        load3 = squeeze(mean(mean(gndTF{1,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load3 = [load3;tmp];

        X(:,1) = [load1;load2;load3];
        X(:,2) = sort(repmat([1 2 3]',subN,1));
        X(:,3) = repmat([1:subN]',3,1);

        [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
    end

    subplot(3,length(timeROI.bins),2*length(timeROI.bins)+t);hold all
    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
    caxis([0, 10])
    h = colorbar;
    %             set(h,'ytick',[0.33:0.03:0.4])
    h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(sprintf('Load-dependent clear-out\nAll:%.1f~%.1fs', timeROI.bins{t}(1),timeROI.bins{t}(2)));
end
saveas(gcf,fullfile(Dir.figs,['RespPowTopoLoadDependent_equal_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
%% to check why there's no correlation overall for young
Nplus1 = load(fullfile(Dir.results,['RespBetaPow',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll');

Np1table = subs;
Np1table.np1= Nplus1.subsAll.betaFrontalACCcorr(:,3);% exmining load 4
[~,Np1table.rank]=sort(Np1table.np1);
[~,Np1table.rank]  = ismember(Np1table.np1,unique(Np1table.np1));

Np1table.beha = mean([dat.beha{2};dat.beha{1}],2);
Np1table.var = mean([dat.betaAvgFrontal{2};dat.betaAvgFrontal{1}],2);
myFigBasic
figure('Position',[300 300 1000 600]);
tmpColor = crameri('grayC',50);
for gi = 1:2
    tmpID = subs.group==gi;
    subplot(1,2,gi);hold all;axis square
            scatter(Np1table.var(tmpID),Np1table.beha(tmpID),2,'k');
            lsline

    for i = 1:sum(tmpID)
        tmp = find(tmpID);
        plot(Np1table.var(tmp(i)),Np1table.beha(tmp(i)), '.','MarkerSize',50,'Color',[tmpColor(Np1table.rank(tmp(i)),:) 0.5]);
        text(Np1table.var(tmp(i))+0.01,Np1table.beha(tmp(i))-0.01,num2str(Np1table.rank(tmp(i))),'HorizontalAlignment','left','Color','r');
    end
title(groupStr{gi})
xlabel('Clear out');ylabel('ACC')
    colormap(tmpColor)
    cb = colorbar;
    cb.Label.String = 'N+1 Rho';
        cb.Label.Rotation = -90;
        cb.Label.Position = cb.Label.Position + [1 0 0];
     cb.Ticks = [0 1];
     cb.TickLabels = {'Low','High'};
end
saveas(gcf,fullfile(Dir.figs,['RespNeuPowCorrLabeled_equal_',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%%
figure('Name','Ftest chans','Position',[100 200 800 300]);
threshP = 0.001;
mkSize = 3;
groupLimit = {[0,90],[0 30]};
for gi = 1:2
    clear X
    for c = 1:length(gndTF{1}.label)

        load1 = squeeze(mean(mean(gndTF{gi,1}.indv(:,c,freqID.beta,timeID.Pre),3),4));%subN
        load2 = squeeze(mean(mean(gndTF{gi,2}.indv(:,c,freqID.beta,timeID.Post),3),4));%subN

        X(:,1) = [load1;load2];
        X(:,2) = sort(repmat([1 2]',sum(subs.group==gi),1));
        X(:,3) = repmat([1:sum(subs.group==gi)]',2,1);

        [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
    end
    subplot(1,3,gi);hold all

    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off');
    caxis(groupLimit{gi})
    h = colorbar;
    h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.bins{t}(1),timeROI.bins{t}(2)),['N=' num2str(sum(subs.group==gi))]);

end

clear X
for c = 1:length(gndTF{1}.label)

    load1 = squeeze(mean(mean(gndTF{1,1}.indv(:,c,freqID.beta,timeID.Pre),3),4));%subN
    tmp = squeeze(mean(mean(gndTF{2,1}.indv(:,c,freqID.beta,timeID.Pre),3),4));%subN
    load1 = [load1;tmp];

    load2 = squeeze(mean(mean(gndTF{1,2}.indv(:,c,freqID.beta,timeID.Post),3),4));%subN
    tmp = squeeze(mean(mean(gndTF{2,2}.indv(:,c,freqID.beta,timeID.Post),3),4));%subN
    load2 = [load2;tmp];

    X(:,1) = [load1;load2];
    X(:,2) = sort(repmat([1 2]',subN,1));
    X(:,3) = repmat([1:subN]',2,1);

    [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
end

subplot(1,3,3);hold all
topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
caxis([0, 25])
h = colorbar;
%             set(h,'ytick',[0.33:0.03:0.4])
h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
h.Label.Rotation = -90;
h.Label.Position = h.Label.Position + [1 0 0];
title(sprintf('All:%.1f~%.1fs', timeROI.bins{t}(1),timeROI.bins{t}(2)));

saveas(gcf,fullfile(Dir.figs,['RespPowTopoTtest_equal_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%% split young and old based on their acc

beha.accAvg = mean([beha.s1_acc,beha.s2_acc, beha.s3_acc],2);
figure('position',[200 200 600 600]);
splitStr = {'HighACC','LowACC'};
threshP = .01;
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
            load1 = squeeze(mean(mean(gndTF{gi,1}.indv(tmp_upper,c,freqID.beta,timeID.Pre),3),4));%subN
            load2 = squeeze(mean(mean(gndTF{gi,2}.indv(tmp_upper,c,freqID.beta,timeID.Post),3),4));%subN
            clear X
            X(:,1) = [load1;load2];
            X(:,2) = sort(repmat([1 2]',sum(tmp_upper),1));
            X(:,3) = repmat([1:sum(tmp_upper)]',2,1);

            [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
        end

        topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off');
        caxis(groupLimit{gi})
        h = colorbar;
        h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.bins{t}(1),timeROI.bins{t}(2)),[splitStr{b},' N=' num2str(sum(tmp_upper))]);
    end
end
saveas(gcf,fullfile(Dir.figs,['RespPowTopoTtestSplitACC_equal_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
