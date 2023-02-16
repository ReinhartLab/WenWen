clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))
addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))

IsCorretTrials = 1;
IsBL2preDelay = 0;
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Young-Old'};
condStr = {'Load 1','Load 2','Load 4'};

frontalROI = {'Fz','F1','F2','AFz','FCz'};
occipROI = {'POz','Oz'};

freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 12];% Hz
timeROI.all = [-5.2 0];% in s
SStime = {[-1.6 timeROI.all(2)],[-2.8 timeROI.all(2)],[-5.2 timeROI.all(2)]};% epoch length for different SS
timeROI.delay = [0 3];% in s
timeROI.pic = {[-1.2 0],[-2.4 -1.2;-1.2 0],[-4.8 -3.6;-3.6 -2.4;-2.4 -1.2;-1.2 0]};
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

%%
freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID.all = dsearchn(gndTF{1}.time',timeROI.all');
timeID.all = timeID.all(1):timeID.all(2);

timeID.delay = dsearchn(gndTF{1}.time',timeROI.delay');
timeID.delay = timeID.delay(1):timeID.delay(2);

for cond_i = 1:3
    for s = 1:size(timeROI.pic{cond_i},1)
        timeID.pic{cond_i,s} = dsearchn(gndTF{1}.time',timeROI.pic{cond_i}(s,:)');
        timeID.pic{cond_i,s} = timeID.pic{cond_i,s}(1):timeID.pic{cond_i,s}(end);
    end
end

clear chanID
[log1,chanID.frontal] = ismember(frontalROI,gndTF{1}.label);
[log2,chanID.occip] = ismember(occipROI,gndTF{1}.label);

clear dat SStimeID

for gi = 1:2
    for cond_i = 1:3
        SStimeID{cond_i} = dsearchn(gndTF{1}.time',SStime{cond_i}');
        SStimeID{cond_i} = SStimeID{cond_i}(1):SStimeID{cond_i}(2);

        dat.betaCurveOccip{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,SStimeID{cond_i}),2),3));
        dat.betaCurveFrontal{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,SStimeID{cond_i}),2),3));

        dat.alphaCurveOccip{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,SStimeID{cond_i}),2),3));
        dat.alphaCurveFrontal{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,SStimeID{cond_i}),2),3));

        for s = 1:size(timeROI.pic{cond_i},1)
            dat.betaAvgFrontal{gi}{cond_i}(:,s) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.pic{cond_i,s}),2),3),4));
            dat.betaAvgOccip{gi}{cond_i}(:,s) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID.pic{cond_i,s}),2),3),4));

            %             dat.alphaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID.delay),2),3),4));
            %             dat.alphaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.delay),2),3),4));
        end
    end
end

% mixed ANOVA: group*cond = 2*3
clear Xa
cond_i  = 3;% examine load 4,group*load
Xa(:,1) = [reshape(dat.betaAvgFrontal{1}{cond_i},size(dat.betaAvgFrontal{1}{cond_i},1)*size(dat.betaAvgFrontal{1}{cond_i},2),1);...
    reshape(dat.betaAvgFrontal{2}{cond_i},size(dat.betaAvgFrontal{2}{cond_i},1)*size(dat.betaAvgFrontal{2}{cond_i},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*4,1);ones(sum(subs.group==2)*4,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;ones(sum(subs.group==1),1)*4;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3;ones(sum(subs.group==2),1)*4];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.betaFrontal, DFs.betaFrontal, MSQs.betaFrontal, Fs.betaFrontal, Ps.betaFrontal]=mixed_between_within_anova(Xa);

T = subs(:,{'name','group','groupStr'});
cond_i  = 3;% examine load 4,group*load
T.s1 = [dat.betaAvgFrontal{2}{cond_i}(:,1);dat.betaAvgFrontal{1}{cond_i}(:,1)];
T.s2 = [dat.betaAvgFrontal{2}{cond_i}(:,2);dat.betaAvgFrontal{1}{cond_i}(:,2)];
T.s3 = [dat.betaAvgFrontal{2}{cond_i}(:,3);dat.betaAvgFrontal{1}{cond_i}(:,3)];
T.s4 = [dat.betaAvgFrontal{2}{cond_i}(:,4);dat.betaAvgFrontal{1}{cond_i}(:,4)];
writetable(T,fullfile(Dir.ana,'TableOutput',['EncodingBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

%% linear fit of SS
clear line_fit
cond_i = 3;
for gi = 1:2
    for sub_i = 1:size(dat.betaAvgFrontal{gi}{cond_i},1)
        line_fit{gi}.params(sub_i,:) = polyfit([1 2 3 4]',dat.betaAvgFrontal{gi}{cond_i}(sub_i,:)',1);
    end
    [~,line_fit{gi}.p,~,line_fit{gi}.Tstat] = ttest(line_fit{gi}.params(:,1));
end

T = subs(:,{'name','group','groupStr'});
T.slope_delay = [line_fit{2}.params(:,1);line_fit{1}.params(:,1)];
% writetable(T,fullfile(Dir.ana,'TableOutput',['EncodingBetaVariance_Slope_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

%% bar & time series
addpath(genpath('D:\Toolbox\crameri_v1.08'))
fig = figure('Position',[100 100 800 400]);
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);

myColors = crameri('batlowW',5);

for gi = 1:2
    subplot(2,2,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveFrontal{gi}{cond_i}));
        se = squeeze(std(dat.betaCurveFrontal{gi}{cond_i},0,1)./sqrt(sum(subs.group==gi)));

        mn2plot = nan(1,length(SStimeID{3}));
        mn2plot(1:length(SStimeID{cond_i})) = mn;
        se2plot = nan(1,length(SStimeID{3}));
        se2plot(1:length(SStimeID{cond_i})) = se;
        shadedErrorBar(gndTF{1}.time(SStimeID{3}),mn2plot,se2plot,{'color',myColors(cond_i,:)})
    end

    legend(condStr,'Location','northeastoutside')
    ytickformat('%.1f')
    ylabel('Variance(AU)');
    xlabel('Time (0s=maintenance)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[frontalROI{:}],freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XTick',[-4.8 -3.6 -2.4 -1.2 0],'YLim',[-1.2 0.3],'XLim',[-5.23 timeROI.all(2)])

    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([-4.8 -4.8],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-3.6 -3.6],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-2.4 -2.4],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([-1.2 -1.2],get(gca,'YLim'),'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k--','HandleVisibility','off')
end

subplot(2,2,3);%bar
hold all;axis square
lineS = {'--','-'};
ms = 10;
for gi = 1:2

    mn = cellfun(@mean,dat.betaAvgFrontal{gi},'UniformOutput',false);
    se = cellfun(@(x)std(x)/sqrt(sum(subs.group==gi)),dat.betaAvgFrontal{gi},'UniformOutput',false);

    errorbar(1:length(mn{1}),mn{1},se{1},'k.','LineWidth',1.5,'LineStyle',lineS{gi},'MarkerSize',ms)

    for c = 1:3
        cond_i = 4-c;
        errorbar(1:length(mn{cond_i}),mn{cond_i},se{cond_i},'.','Color',myColors(cond_i,:),'LineWidth',1.5,'LineStyle',lineS{gi},'MarkerSize',ms,'HandleVisibility','off')
        % errorbar(1:length(mn{cond_i}),mn{cond_i},se{cond_i},'.','Color',groupColor(gi,:),'LineWidth',1.5,'LineStyle',lineS{gi},'MarkerSize',ms)
    end
end
legend(groupStr,'Location','bestoutside')
xlabel('Set size')
set(gca,'XLim',[0 5],'XTick',1:4)
ylabel('Variance(AU)')
title('beta encoding')

% plot slope index

subplot(2,2,4);hold all;box on;axis square
subs.slope = [line_fit{2}.params(:,1);line_fit{1}.params(:,1)];
boxplot(subs.slope,subs.groupStr)
plot(get(gca,'xlim'),[0 0],'k--')
set(gca,'YLim', [-0.4 0.2],'ytick',-0.4:0.2:0.6,'XLim',[0 3]);
ylabel({'\bfSlope'});

for gi = 1:2
    if  line_fit{gi}.p<.05
        text(gi,0.1,'*','HorizontalAlignment','center','FontSize',20)
    end
end
[~,p] = ttest2(line_fit{1}.params(:,1),line_fit{2}.params(:,1));
text(1.5,-0.35,['p=',num2str(p)],'HorizontalAlignment','center','FontSize',12);
plot([1 2],[-0.3 -0.3],'k')

saveas(gca,fullfile(Dir.figs,['EncodingBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);

behaStr = {'RT','Accuracy'};

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',5);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

var_Str = {'average','last item'};

for var_idx = 1:2
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
            for cond = 1:3
                if var_idx ==1
                    X = mean(dat.betaAvgFrontal{gi}{cond},2);
                else
                    X = dat.betaAvgFrontal{gi}{cond}(:,end);
                end
                [r,pval] = corr(X,dat.beha{gi}(:,cond),'type','Spearman');

                Y = dat.beha{gi}(:,cond);

                g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
                g(gi,cond).geom_point(); %Draw a scatter plot
                g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots

                g(gi,cond).set_names('x',['Beta variance of ' var_Str{var_idx}],'y',behaStr{idx}); %Set the title of the axis
                g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
                g(gi,cond).set_title(sprintf('%s: Rho = %.3f, p = %.3f',groupStr{gi},r,pval));
                g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
            end

        end
        fig = figure('Position',[100 100 800 500]);
        g.draw();
        saveas(gcf,fullfile(Dir.figs,['EncodingNeuVarCorr_',behaStr{idx},var_Str{var_idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
    end
end