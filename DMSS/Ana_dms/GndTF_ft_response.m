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
IsBL2preDelay = 1;
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
groupStr = {'Young','Old','Young-Old'};
condStr = {'ss1','ss2','ss4'};
condDiffStr = {'S2-S1','S4-S1','S4-S2'};

frontalROI = {'Fz','F1','F2','FCz','FC1','FC2','C1','C2','Cz'};
frontalROI = {'Fz','F1','F2'};
occipROI = {'Pz','POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz
timeROI = [-0.4 0.8];% in s, for curve plotting
timeROIanova = [0 0.8];% in s
%%
for sub_i = 1:subN
    subname = subs.name{sub_i};

    matName = fullfile(Dir.results,[subname,'_resp_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
    if isfile(matName)
        load(matName);
        for s = 1:3
            subsAll{sub_i,s} = ft_freqdescriptives([],tfDat{s});
        end
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
        cfg.xlim  = [-0.4 0.8];
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
        rectangle('Position',[0 0.5 0.8 40],'LineWidth',1.2,'LineStyle','--')

        subplot(2,3,3+cond_i);axis square

        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROIanova;%time
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

        ft_topoplotTFR(cfg, gndTF{gi,cond_i});
        colorbar

        % %         ft_singleplotTFR(cfg, gndTFdiff{gi,cond_i});
        % %         rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')
        % %         title(condDiffStr{cond_i},groupStr{gi},'FontSize',14);
        % %
        % %         hc = colorbar;
        % %         hc.Label.String = 'Power(dB)';
        % %         hc.Label.Rotation = -90;
        % %         hc.Label.Position = [4.2 0 0];

    end
    saveas(gcf,fullfile(Dir.figs,['Resp_power_',groupStr{gi},[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4}, '.bmp']))
end

%%  time frequency of group average
cfg = [];
cfg.xlim  = [-0.4 0.8];
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

figure('Position',[100 100 800 290]);
for gi = 1:3
    subplot(1,3,gi);    axis square;hold all

    ft_singleplotTFR(cfg, GroupAvg{gi});

    title(groupStr{gi},'Average set sizes','FontSize',14);
    hc = colorbar;
    hc.Label.String = 'Power(dB)';
    hc.Label.Rotation = -90;
    hc.Label.Position = [4.2 0 0];

    xlabel('Time(0s=response)')
    rectangle('Position',[0 -0.5 0.8 50],'LineWidth',1.2,'LineStyle','--')
    set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',timeROI,'YLim',[0.5 40])

end
saveas(gcf,fullfile(Dir.figs,['Resp_power_GroupAvg_',[cfg.channel{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4}, '.bmp']))

%% topography

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

for gi = 1:2
    figure('Position',[100 100 900 400])
    for cond_i = 1:3
        subplot(2,3,cond_i);axis square
        ft_topoplotTFR(cfg, gndTF{gi,cond_i});
        title(condStr{cond_i},groupStr{gi});colorbar

        subplot(2,3,cond_i+3);

        ft_topoplotTFR(cfg, gndTFdiff{gi,cond_i});
        title(condDiffStr{cond_i},groupStr{gi});colorbar

    end
    %     saveas(gcf,fullfile(Dir.figs,['TimeFreq_resp',groupStr{gi},num2str(cfg.xlim(1)),'~',num2str(cfg.xlim(2)),'s',[cfg.highlightchannel{:}] '_topo.bmp']))
end

%%

freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID = dsearchn(gndTF{1}.time',timeROI');
timeID = timeID(1):timeID(2);

timeIDanova = dsearchn(gndTF{1}.time',timeROIanova');
timeIDanova = timeIDanova(1):timeIDanova(2);

clear chanID
[log1,chanID.frontal] = ismember(frontalROI,gndTF{1}.label);
[log2,chanID.occip] = ismember(occipROI,gndTF{1}.label);

clear dat

for gi = 1:2
    for cond_i = 1:3
        dat.betaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeIDanova),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.betaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeIDanova),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        %         dat.alphaAvg{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeIDanova),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        dat.betaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID),2,"omitnan"),3,"omitnan"));
        dat.betaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID),2,"omitnan"),3,"omitnan"));

        %         dat.alphaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID),2,"omitnan"),3,"omitnan"));
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


%% bar plot
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
title(sprintf('Frontal beta(%.1f~%.1fs)',timeROIanova(1),timeROIanova(2)))
% plot significance
if Ps.betaFrontal{1}<.05
    plot([1 2],[1.2 1.2],'k','HandleVisibility','off')
    sigAster = '*';
    if Ps.betaFrontal{1}<.01
        sigAster = '**';
    elseif Ps.betaFrontal{1}<.001
        sigAster = '***';
    end
    text(1.5,1.3,sigAster,'FontSize',18,'HorizontalAlignment','center')
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
title(sprintf('Posterior beta(%.1f~%.1fs)',timeROIanova(1),timeROIanova(2)))

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
    legend(condStr,'Location','southeast')
    ytickformat('%.1f')
    set(gca,'xlim',timeROI)
    ylabel('Power(dB)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},'Frontal beta')

    subplot(3,2,gi+2);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    legend(condStr)
    ytickformat('%.1f')
    ylabel('Power(dB)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},'Occipital beta')
end

saveas(gca,fullfile(Dir.figs,['RespBetaPower_',num2str(timeROIanova(1)),'~',num2str(timeROIanova(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);
dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));

% dat.beha{1} = table2array(beha(beha.group==1,{'s1_d','s2_d','s3_d'}));
% dat.beha{2} = table2array(beha(beha.group==2,{'s1_d','s2_d','s3_d'}));
%
% dat.beha{1} = table2array(beha(beha.group==1,{'s1_rt','s2_rt','s3_rt'}));
% dat.beha{2} = table2array(beha(beha.group==2,{'s1_rt','s2_rt','s3_rt'}));

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

clear g
for gi = 1:2
    [r,pval] = corr(dat.betaAvgOccip{gi},dat.beha{gi},'type','Spearman');

    for cond = 1:3
        Y = dat.beha{gi}(:,cond);
        X = dat.betaAvgOccip{gi}(:,cond);

        g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,cond).geom_point(); %Draw a scatter plot
        g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
        g(gi,cond).set_names('x','Beha power(dB)','y','Beha'); %Set the title of the axis
        g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
        g(gi,cond).set_title(sprintf('%s: Rho = %.3f, p = %.3f',groupStr{gi},r(cond,cond),pval(cond,cond)));
        g(gi,cond).set_text_options('base_size' ,8,'title_scaling' ,1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    end
end
figure('Position',[100 100 1000 600]);
g.draw();

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
    xlabel('Time(0s=response)')
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
    xlabel('Time(0s=response)')
    title(groupStr{gi},'Posterior beta')

end
saveas(gca,fullfile(Dir.figs,['RespBetaPowerBehavGLM_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.bmp']))
