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
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
groupStr = {'Young','Old','Young-Old'};
condStr = {'ss1','ss2','ss4'};
condDiffStr = {'S2-S1','S4-S1','S4-S2'};

frontalROI = {'Fz','F1','F2','FCz','FC1','FC2'};
% frontalROI = {'Fz','F1','F2',};
% frontalROI = {'Fz','F1','F2','FCz','FC1','FC2','Cz','C1','C2'};
occipROI = {'Pz','POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz
timeROI = [-1 3];% in s
%%
for sub_i = 1:subN
    subname = subs.name{sub_i};

    matName = fullfile(Dir.results,[subname,'_delay_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
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
    figure('Position',[100 100 800 400]);
    axis square

    for cond_i = 1:3
        subplot(2,3,cond_i);

        cfg = [];
        cfg.xlim  = [-1 4.5];
        cfg.ylim = [1 40];
        cfg.zlim = [-1.5 1.5];
        % cfg.colormap = bkr;
        cfg.colormap = jet;
        cfg.channel = frontalROI;
        cfg.figure = 'gca';
        % cfg.maskparameter = 'mask';
        % cfg.maskstyle = 'outline';
        % cfg.maskstyle = 'opacity';
        % cfg.maskalpha = 0.3;

        ft_singleplotTFR(cfg, gndTF{gi,cond_i});

        title(condStr{cond_i},groupStr{gi},'FontSize',14);
        hc = colorbar;
        hc.Label.String = 'Power(dB)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',timeROI)
        if cond_i == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time(0s=maintenance)')
        rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')

        subplot(2,3,3+cond_i);axis square

        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROI;%time
        cfg.zlim = [-1 1];
        cfg.highlightchannel = frontalROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';

        ft_topoplotTFR(cfg, gndTF{gi,cond_i});
        colorbar

    end
    saveas(gcf,fullfile(Dir.figs,['TimeFreq_',groupStr{gi},[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))
end

%%  time frequency of group average
cfg = [];
cfg.xlim  = [-1 5];
cfg.ylim = [1 40];
cfg.zlim = [-1.5 1.5];
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
    subplot(1,3,gi);    axis square;hold all

    ft_singleplotTFR(cfg, GroupAvg{gi});

    title(groupStr{gi},'Average SS','FontSize',14);
    hc = colorbar;
    hc.Label.String = 'Power(dB)';
    hc.Label.Rotation = -90;
    hc.Label.Position = [4.2 0 0];

    xlabel('Time(0s=maintenance)')
    rectangle('Position',[0 -0.4 3 50],'LineWidth',1.2,'LineStyle','--')
    set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',timeROI,'YLim',[0.5 40])
end
saveas(gcf,fullfile(Dir.figs,['TimeFreq_GroupAvg_',[cfg.channel{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))

%% time frequency of group diff

figure('Position',[100 100 800 400])
for cond_i = 1:3
    subplot(2,3,cond_i);axis square
    cfg = [];
    cfg.xlim  = [-1 4.5];
    cfg.ylim = [1 40];
    cfg.zlim = [-1.5 1.5];
    % cfg.colormap = bkr;
    cfg.colormap = jet;

    % cfg.channel = {'all'};
    cfg.channel = frontalROI;
    cfg.figure = 'gca';
    % cfg.maskparameter = 'mask';
    % cfg.maskstyle = 'outline';
    % cfg.maskstyle = 'opacity';
    % cfg.maskalpha = 0.3;

    ft_singleplotTFR(cfg, gndTF_Groupdiff{cond_i});
    rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')
    title(condStr{cond_i},groupStr{3},'FontSize',14);colorbar

    subplot(2,3,cond_i+3);
    cfg = [];
    cfg.ylim = freq.betaFreq;%freq
    cfg.xlim = timeROI;%time
    cfg.zlim = [-1 1];
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
saveas(gcf,fullfile(Dir.figs,['TimeFreqGroupDiff_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))

%%

freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID = dsearchn(gndTF{1}.time',timeROI');
timeID = timeID(1):timeID(2);

clear chanID
[log1,chanID.frontal] = ismember(frontalROI,gndTF{1}.label);
[log2,chanID.occip] = ismember(occipROI,gndTF{1}.label);

clear dat

for gi = 1:2
    for cond_i = 1:3
        dat.betaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID),2),3),4));
        dat.betaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID),2),3),4));

        dat.alphaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID),2),3),4));
        dat.alphaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID),2),3),4));

        dat.betaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID),2),3));
        dat.betaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID),2),3));

        dat.alphaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID),2),3));
        dat.alphaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID),2),3));
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


%% bar plot -beta
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[100 100 700 1000]);

subplot(3,2,5);hold all;axis square
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
errorbar(xcord,mn,se,'k.')

legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Power(dB)')
title(sprintf('Frontal beta(%.1f~%.1fs)',timeROI(1),timeROI(2)))

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

subplot(3,2,6);hold all;axis square
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
ylabel('Power(dB)')
title(sprintf('Posterior beta(%.1f~%.1fs)',timeROI(1),timeROI(2)))

myColors = crameri('batlowW',5);

times = gndTF{1}.time(timeID);
for gi = 1:2
    subplot(3,2,gi);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveFrontal{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveFrontal{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Power(dB)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Frontal %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))

    subplot(3,2,gi+2);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Power(dB)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
end

saveas(gca,fullfile(Dir.figs,['DelayBetaPower_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))
%%
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[100 100 700 1000]);

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
ylabel('Power(dB)')
title(sprintf('Posterior alpha(%.1f~%.1fs)',timeROI(1),timeROI(2)))

subplot(3,2,6);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';
% hb = bar(mn,'FaceAlpha',0.8);
% xcord = vertcat(hb(:).XEndPoints)';

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
ylabel('Power(dB)')

title(sprintf('Frontal alpha(%.1f~%.1fs)',timeROI(1),timeROI(2)))
% plot significance
% plot(xcord(1,[2 3]),[1.2 1.2],'k','HandleVisibility','off')
% text(mean(xcord(1,[2 3])),1.2,'***','FontSize',18,'HorizontalAlignment','center')
% plot(xcord(1,[1 3]),[1 1],'k','HandleVisibility','off')
% text(mean(xcord(1,[1 3])),1,'***','FontSize',18,'HorizontalAlignment','center')

myColors = crameri('batlowW',5);

times = gndTF{1}.time(timeID);
for gi = 1:2
    subplot(3,2,gi);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.alphaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Power(dB)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.alphaFreq(1),freq.alphaFreq(2)))

    subplot(3,2,gi+2);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveFrontal{gi}(:,cond_i,:)));
        se = squeeze(std(dat.alphaCurveFrontal{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Power(dB)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Frontal %d~%dHz',freq.alphaFreq(1),freq.alphaFreq(2)))
end

saveas(gca,fullfile(Dir.figs,['DelayAlphaPower_',num2str(freq.alphaFreq(1)),'~',num2str(freq.alphaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);
% dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
% dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));

dat.beha{1} = table2array(beha(beha.group==1,{'s1_d','s2_d','s3_d'}));
dat.beha{2} = table2array(beha(beha.group==2,{'s1_d','s2_d','s3_d'}));

% dat.beha{1} = table2array(beha(beha.group==1,{'s1_rt','s2_rt','s3_rt'}));
% dat.beha{2} = table2array(beha(beha.group==2,{'s1_rt','s2_rt','s3_rt'}));

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',4);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

clear g
for gi = 1:2
    [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Spearman');

    for cond = 1:3
        Y = dat.beha{gi}(:,cond);
        X = dat.betaAvgFrontal{gi}(:,cond);

        g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,cond).geom_point(); %Draw a scatter plot
        g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
        g(gi,cond).set_names('x',sprintf('Beta power(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2)),'y','Beha'); %Set the title of the axis
        g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
        g(gi,cond).set_title(sprintf('%s SS#%d\nRho = %.3f, p = %.3f',groupStr{gi},cond,r(cond,cond),pval(cond,cond)));
        g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    end

    cond = 4;
    Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
    X = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);
    [r,pval] = corr(Y,X,'type','Spearman');

    g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
    g(gi,cond).geom_point(); %Draw a scatter plot
    g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
    g(gi,cond).set_names('x',sprintf('Beta power(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2)),'y','Beha'); %Set the title of the axis
    g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
    g(gi,cond).set_title(sprintf('%s collapse\nRho = %.3f, p = %.3f',groupStr{gi},r,pval));
    g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times


end
figure('Position',[100 100 1100 550]);
g.draw();
saveas(gcf,fullfile(Dir.figs,['DelayPowerCorr_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))

%% GLM beta coeff
figure('Position',[500 500 650 600]);
threshP = .05;
for gi = 1:2
    subplot(2,2,gi);hold all;axis square;
    for cond_i = 1:3
        X = dat.beha{gi}(:,cond_i);
        for t = 1:size(dat.betaCurveFrontal{gi},3)
            Y = squeeze(dat.betaCurveFrontal{gi}(:,cond_i,t));
            [b,dev,stats] = glmfit(X,Y,'normal');
            glmBeta.b(gi,cond_i,t) = b(2);
            glmBeta.p(gi,cond_i,t) = stats.p(2);
        end

        plot(times,squeeze(glmBeta.b(gi,cond_i,:)),'Color',myColors(cond_i,:))

        mySig = nan(length(times),1);
        mySig(glmBeta.p(gi,cond_i,:)<threshP) = min(glmBeta.b(gi,cond_i,:))*1.1;
        plot(times,mySig,'Color',myColors(cond_i,:),'LineWidth',3,'HandleVisibility','off')
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('\beta coefficient(power~accuracy)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},'Frontal beta')

    subplot(2,2,gi+2);hold all;axis square;

    for cond_i = 1:3
        X = dat.beha{gi}(:,cond_i);
        for t = 1:size(dat.betaCurveOccip{gi},3)
            Y = squeeze(dat.betaCurveOccip{gi}(:,cond_i,t));
            [b,dev,stats] = glmfit(X,Y,'normal');
            glmBeta.b(gi,cond_i,t) = b(2);
            glmBeta.p(gi,cond_i,t) = stats.p(2);
        end

        plot(times,squeeze(glmBeta.b(gi,cond_i,:)),'Color',myColors(cond_i,:))

        mySig = nan(length(times),1);
        mySig(glmBeta.p(gi,cond_i,:)<threshP) = min(glmBeta.b(gi,cond_i,:))*1.1;
        plot(times,mySig,'Color',myColors(cond_i,:),'LineWidth',3,'HandleVisibility','off')
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('\beta coefficient(power~accuracy)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},'Posterior beta')

end
saveas(gca,fullfile(Dir.figs,['DelayBetaPowerBehavGLM_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))
% % % % %%  cluster for each condition and each group
% % % % % if exist(permFile,'file') == 0
% % % % clear perm
% % % % for gi = 1:2
% % % %     subN = length(mySubs{gi});
% % % %     for cond = 1:3
% % % %         clear datP datDiffP
% % % %
% % % %         datP = gndTF_SSdiff{gi,cond}.powspctrm(:,chanID.beta,:,:);% subN*chan*freq*time
% % % %         datDiffP = gndTF_SSdiff{gi,cond}.powspctrm(:,chanID.beta,:,:);
% % % %
% % % %         datP = squeeze(mean(datP,2));
% % % %         datP = permute(datP,[2 3 1]);% freq*time*subs
% % % %
% % % %         datDiffP = squeeze(mean(datDiffP,2));
% % % %         datDiffP = permute(datDiffP,[2 3 1]);% freq*time*subs
% % % %
% % % %         %             BLidx = dsearchn(tfAll(1).all{1}.time',[-0.4 0]');
% % % %         %             datP = bsxfun(@rdivide,datP,mean(datP(:,BLidx(1):BLidx(2),:),2))-1;
% % % %
% % % %         nperm = 10000;
% % % %
% % % %         diffstat = 'diff';%'t';% or 'diff'
% % % %         [datobs, datrnd] = cluster_test_helper(datP, nperm, diffstat);
% % % %         tail = 0;alpha = .05;
% % % %         clusteralpha = .05;
% % % %         clusterstat = 'sum'; %'size'
% % % %
% % % %         [perm(gi).pval(cond,:,:), perm(gi).h(cond,:,:), perm(gi).clusterinfo{cond}] = cluster_test(datobs, datrnd, tail, alpha,...
% % % %             clusteralpha, clusterstat);
% % % %
% % % %         [datobs, datrnd] = cluster_test_helper(datDiffP, nperm, diffstat);
% % % %         [perm(gi).pvalDiff(cond,:,:), perm(gi).hDiff(cond,:,:), perm(gi).clusterinfoDiff{cond}] = cluster_test(datobs, datrnd, tail, alpha,...
% % % %             clusteralpha, clusterstat);
% % % %     end
% % % % end
% % % % save(permFile,'perm')
% % % % % else
% % % % %     load(permFile)
% % % % % end
% % % % %% frontal mask
% % % % w = [1 1 1];
% % % % b = [0,0,1];
% % % % k = [0,0,0];
% % % % r = [1,0,0];
% % % % bkr = createcolormap(w,b,k,r,w); % 256x3 array
% % % %
% % % % cfg = [];
% % % % cfg.xlim  = [-1 3];
% % % % cfg.ylim = [1 40];
% % % % cfg.zlim = [-1.5 1.5];
% % % % % cfg.colormap = bkr;
% % % % cfg.colormap = jet;
% % % %
% % % % cfg.channel = frontalROI;
% % % % cfg.figure = 'gca';
% % % % % cfg.maskparameter = 'mask';
% % % % % cfg.maskstyle = 'outline';
% % % %
% % % % % cfg.maskstyle = 'opacity';
% % % % % cfg.maskalpha = 0.3;
% % % %
% % % % for gi = 1:2
% % % %     figure('Position',[100 100 800 400]);
% % % %     axis square
% % % %
% % % %     for cond_i = 1:3
% % % %         subplot(2,3,cond_i);
% % % %
% % % %         tmpmask = squeeze(perm(gi).h(cond_i,:,:));%create mask based on permutation
% % % %         tmpmask = repmat(tmpmask,1,1,length(gndTF{gi,cond_i}.label));
% % % %         gndTF{gi,cond_i}.mask = permute(tmpmask,[3 1 2]);
% % % %
% % % %         ft_singleplotTFR(cfg, gndTF{gi,cond_i});
% % % %
% % % %         title(condStr{cond_i},groupStr{gi},'FontSize',14);
% % % %         hc = colorbar;
% % % %         hc.Label.String = 'Power(dB)';
% % % %         hc.Label.Rotation = -90;
% % % %         hc.Label.Position = [4.2 0 0];
% % % %
% % % %         set(gca,'ytick',0:10:length(gndTF{1}.freq))
% % % %         if cond_i == 1
% % % %             ylabel('\bfFrequency(Hz)');
% % % %         end
% % % %         %         xlabel('Time(0s=maintenance)')
% % % %         rectangle('Position',[0 0.5 2 40],'LineWidth',1.2,'LineStyle','--')
% % % %
% % % %         subplot(2,3,3+cond_i);
% % % %         tmpmask = squeeze(perm(gi).hDiff(cond_i,:,:));%create mask based on permutation
% % % %         tmpmask = repmat(tmpmask,1,1,length(gndTF_SSdiff{gi,cond_i}.label));
% % % %         gndTF_SSdiff{gi,cond_i}.mask = permute(tmpmask,[3 1 2]);
% % % %
% % % %         ft_singleplotTFR(cfg, gndTF_SSdiff{gi,cond_i});
% % % %         rectangle('Position',[0 0.5 2 40],'LineWidth',1.2,'LineStyle','--')
% % % %         title(condDiffStr{cond_i},groupStr{gi},'FontSize',14);
% % % %
% % % %         hc = colorbar;
% % % %         hc.Label.String = 'Power(dB)';
% % % %         hc.Label.Rotation = -90;
% % % %         hc.Label.Position = [4.2 0 0];
% % % %
% % % %         set(gca,'ytick',0:10:length(gndTF{1}.freq))
% % % %         if cond_i == 1
% % % %             ylabel('\bfFrequency(Hz)');
% % % %         end
% % % %         xlabel('Time(0s=maintenance)')
% % % %     end
% % % %     %     saveas(gcf,fullfile('Figs',['Dots_',groupStr{gi},[cfg.channel{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.bmp']))
% % % % end
