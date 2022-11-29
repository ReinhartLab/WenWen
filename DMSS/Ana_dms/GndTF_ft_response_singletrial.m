clear;
rng('shuffle');
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

        dat.alphaAvg{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeIDanova),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        dat.betaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID),2,"omitnan"),3,"omitnan"));
        dat.betaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID),2,"omitnan"),3,"omitnan"));

        dat.alphaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID),2,"omitnan"),3,"omitnan"));
    end
end


%% bar plot
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[100 100 700 1000]);

subplot(3,2,5);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==1))]';

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
plot(xcord(1,:),dat.alphaAvg{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvg{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

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
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==1))]';

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

saveas(gca,fullfile(Dir.figs,['RespBetaPower_',num2str(timeROIanova(1)),'~',num2str(timeROIanova(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))

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
saveas(gca,fullfile(Dir.figs,['RespBetaPowerBehavGLM_',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
