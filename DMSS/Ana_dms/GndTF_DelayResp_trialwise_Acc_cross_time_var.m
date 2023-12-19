clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))
addpath(genpath('D:\intWM-E\toolbox\distributionPlot'))

IsCorretTrials = 0;% enforced for accuracy
IsLap = 1;
IsdePhase=1;
IsBL2preDelay = 1;% enforced

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Older','Both'};
condStr = {'ss1','ss2','ss4'};
ssN = [1 2 4];

% chanROI = {'CPz','CP1','CP2','Pz','Cz'};
chanROI = {'FCz','Fz','F1','F2','AFz'};

freq.betaFreq = [15 25];% Hz
% freq.betaFreq = [8 12];% Hz
% freq.betaFreq = [1 3];% Hz
% freq.betaFreq = [4 7];% Hz
% freq.betaFreq = [28 40];% Hz

timeROI.Delay = [0 3];% in s
% timeROI.Post = [0 0.6];% in s
timeROI.Post = [0.1 0.5];% in s
timeROI.bl = [-0.4 -0.1]; 
%%
tbl  = table; % for maintenance analysis,table contains RT(n), delay power(n)
tblN_1 = table;% for n-1 analysis, table contains RT(n), delay power(n), deletion power(n-1)

for sub_i = 1:subN
    subname = subs.name{sub_i};
    disp(subname)
    matName =fullfile(Dir.results,[subname,'_delayresp_trialwise',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);

        freqID.beta = dsearchn(tfDat.PowN.freq',freq.betaFreq');
        freqID.beta = freqID.beta(1):freqID.beta(2);

        timeID.Post = dsearchn(tfDat.PowN_1.time',timeROI.Post');
        timeID.Post = timeID.Post(1):timeID.Post(2);

        timeID.bl4Post = dsearchn(tfDat.PowN_1.time',timeROI.bl');
        timeID.bl4Post = timeID.bl4Post(1):timeID.bl4Post(2);

        timeID.bl4Delay = dsearchn(tfDat.PowN.time',timeROI.bl');
        timeID.bl4Delay = timeID.bl4Delay(1):timeID.bl4Delay(2);

        timeID.Delay = dsearchn(tfDat.PowN.time',timeROI.Delay');
        timeID.Delay = timeID.Delay(1):timeID.Delay(2);

        clear chanID
        [log1,chanID.frontal] = ismember(chanROI,tfDat.PowN.label);

        M.subj = ones(height(M),1)*sub_i;
        M.group = ones(height(M),1)*subs.group(sub_i);

    % ----cross time variance
    M.betaAvgM = squeeze(mean(mean(var(tfDat.PowN.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.Delay),0,4,'omitnan'),2,"omitnan"),3,"omitnan"));
    M.betaAvgD = squeeze(mean(mean(var(tfDat.PowN_1.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.Post),0,4,'omitnan'),2,"omitnan"),3,"omitnan"));

%     M.betaAvgM_bl = squeeze(mean(mean(var(tfDat.PowN.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.bl4Delay),0,4,'omitnan'),2,"omitnan"),3,"omitnan"));
%     M.betaAvgD_bl = squeeze(mean(mean(var(tfDat.PowN_1.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.bl4Post),0,4,'omitnan'),2,"omitnan"),3,"omitnan"));
% 
%     M.betaAvgM = log(M.betaAvgM./M.betaAvgM_bl);
%     M.betaAvgD = log(M.betaAvgD./M.betaAvgD_bl);

        tbl = [tbl;M];%

        %----------------
        tmp_cidx = ismember(M.trl-1,M.trl);
        Mc = M(tmp_cidx,:);

        tmp_pidx = ismember(M.trl,Mc.trl-1);
        Mp = M(tmp_pidx,{'betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'});
        Mp = renamevars(Mp,{'betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'},{'betaAvgD_1','trl_1','rt_1','corr_1','ss_1'});

        Mc = [Mc Mp];
        tblN_1 = [tblN_1;Mc];    
    end
end

save(fullfile(Dir.results,['DelayRespBetaPow_trialwise_acc_crosstime_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'tbl','tblN_1','-v7.3')

%% glme for frontal beta power 
load(fullfile(Dir.results,['DelayRespBetaPow_trialwise_acc_crosstime_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))
writetable(tblN_1,fullfile(Dir.ana,'TableOutput',['DelayRespBetaPow_trialwise_acc_crosstime_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))
%--------selected channels
tblN_1.subj = categorical(tblN_1.subj);
tblN_1.ss_num = categorical(tblN_1.ss_num); %n ss
tblN_1.group = categorical(tblN_1.group);
tblN_1.ss_1 = categorical(tblN_1.ss_1);% n-1 ss

glme = fitglme(tblN_1,'button_resp_corr ~ betaAvgM*betaAvgD_1*group*ss_num+(1+betaAvgM+betaAvgD_1+ss_num|subj)','Distribution','Binomial','Link','logit');

t = dataset2table(glme.anova);
t.formula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_crosstime_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_anova1'])
t = dataset2table(glme.Coefficients);
t.formula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_crosstime_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_Coefficients1'])

return
%% note: the mean of db normalized single trial data might be different from db normalization on trial-average data 

% a = [-rand(15,1);-200]*10;% single trial power
% bl = -10;% baseline power
% mean(10*log10(a/bl))
% 10*log10(mean(a)/bl)

%%
myFigBasic
myColor = [0 0 0;3 154 53;102 0 204]./255;
lineW = 1.2;
fig = figure('Position',[200 200 1600 300]);

tblstats = grpstats(tblN_1,['button_resp_corr'],["mean","sem"],"DataVars",["betaAvgM","betaAvgD"]);
subplot(1,4,1);hold all;axis square
errorbar([0 1],tblstats.mean_betaAvgM,tblstats.sem_betaAvgM,'k')
set(gca,'XLim',[-1 2],'XTick',[0 1],'XTickLabel',{'Incorrect','Correct'})
xlabel('Memory accuracy of trial N')
ylabel('Maintenance beta power of trials N')
title([chanROI{:}],'n-n')

subplot(1,4,2);hold all;axis square
errorbar([0 1],tblstats.mean_betaAvgD,tblstats.sem_betaAvgD,'k')
set(gca,'XLim',[-1 2],'XTick',[0 1],'XTickLabel',{'Incorrect','Correct'})
xlabel('Memory accuracy of trial N')
ylabel('Deletion beta power of trials N')
title([chanROI{:}],'n-n')

subplot(1,4,3);hold all;axis square
halfStr = {'left','right'};
tblstats = grpstats(tblN_1,["group",'button_resp_corr'],["mean","sem"],"DataVars","betaAvgD_1");
for gi = 1:2
    distributionPlot(tblN_1.betaAvgD_1(double(tblN_1.group)==gi & tblN_1.button_resp_corr ==0),'xValues',0+(gi-1.5)*0.5,'histOri',halfStr{gi},'color',myColor(gi+1,:),'widthDiv',[2 2],'showMM',2)
    distributionPlot(tblN_1.betaAvgD_1(double(tblN_1.group)==gi & tblN_1.button_resp_corr ==1),'xValues',1+(gi-1.5)*0.5,'histOri',halfStr{gi},'color',myColor(gi+1,:),'widthDiv',[2 2],'showMM',2)

    tmp_tbl = tblN_1(double(tblN_1.group)==gi,:);

%     corArray = tmp_tbl.betaAvgD_1(double(tmp_tbl.group)==gi & tmp_tbl.button_resp_corr==1);
%     incorArray = tmp_tbl.betaAvgD_1(double(tmp_tbl.group)==gi & tmp_tbl.button_resp_corr==0);
end
set(gca,'XLim',[-1 2],'XTick',[0 1],'XTickLabel',{'Incorrect','Correct'})
set(gca,'YLim',[-11 20],'YTick',-10:10:20)
xlabel('Memory accuracy of trial N')
ylabel('Deletion beta power of trials N-1')
legend(groupStr,'Location','bestoutside')
title([chanROI{:}],'n+1')

axes('Position',[.57 .7 .05 .15]);hold all
for gi = 1:2
    errorbar([0 1],tblstats.mean_betaAvgD_1(double(tblstats.group)==gi),tblstats.sem_betaAvgD_1(double(tblstats.group)==gi),'Color',myColor(gi+1,:),'LineWidth',lineW)

    tmp_tbl = tblN_1(double(tblN_1.group)==gi,:);
    corArray = tmp_tbl.betaAvgD_1(double(tmp_tbl.group)==gi & tmp_tbl.button_resp_corr==1);
    incorArray = tmp_tbl.betaAvgD_1(double(tmp_tbl.group)==gi & tmp_tbl.button_resp_corr==0);

    [gstats.pval(gi), gstats.t(gi),gstats.df(gi)] = mult_comp_perm_t2(corArray,incorArray,10000,0,0.001);
    gstats.d(gi) = computeCohen_d(corArray,incorArray);
    if gstats.pval(gi)<.05
        plot([0 1],[1 1]*mean(incorArray)*1.3,'Color',myColor(gi+1,:),'HandleVisibility','off')
        text(0.5,mean(incorArray)*1.3,sprintf('*\np = %.3f',gstats.pval(gi)),'Color',myColor(gi+1,:),'HorizontalAlignment','center')
    end
end

set(gca,'XLim',[-0.5 1.5],'XTick',[0 1],'XTickLabel',{'Incorrect','Correct'})
set(gca,'YLim',[-1.2 -0.4],'YTick',[-1.2 -0.4])
ytickformat('%.1f')
%%
subplot(1,4,4);hold all;
tblstats = grpstats(tblN_1,["group",'button_resp_corr','ss_num'],["mean","sem"],"DataVars","betaAvgD_1");

for gi = 1:2
    for s = 1:3
        errorbar([0 1]+gi*2,tblstats.mean_betaAvgD_1(double(tblstats.group)==gi & double(tblstats.ss_num) == s),tblstats.sem_betaAvgD_1(double(tblstats.group)==gi & double(tblstats.ss_num) == s),'Color',myColor(s,:),'LineWidth',lineW)

        tmp_tbl = tblN_1(double(tblN_1.group)==gi & double(tblN_1.ss_num)==s,:);

        corArray = tmp_tbl.betaAvgD_1(double(tmp_tbl.group)==gi & tmp_tbl.button_resp_corr==1);
        incorArray = tmp_tbl.betaAvgD_1(double(tmp_tbl.group)==gi & tmp_tbl.button_resp_corr==0);

%         [tstats.pval(gi,s), tstats.t(gi,s),tstats.df(gi,s)] = mult_comp_perm_t2(corArray,incorArray,10000);
%         tstats.d(gi,s) = computeCohen_d(corArray,incorArray);
%         if tstats.pval(gi,s)<.05
%             plot([0 1]+gi*2,[1 1]*mean(incorArray)*1.3,'Color',myColor(s,:))
%             text(0.5+gi*2,mean(incorArray)*1.3,sprintf('*\np = %.3f',tstats.pval(gi,s)),'Color',myColor(s,:),'HorizontalAlignment','center')
%         end
    end
    text(0.5+gi*2,-2,groupStr{gi},'HorizontalAlignment','center')
end
 set(gca,'XLim',[1 6],'XTick',[2:5],'XTickLabel',{'Incorrect','Correct'})
    xlabel('Memory accuracy of trial N')
ylabel('Deletion beta power of trials N-1')
legend(condStr,'Location','bestoutside')
title([chanROI{:}],'n+1: ss')

% saveas(gca,fullfile(Dir.figs,['DelayRespPowerTrial_glme_ACC_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
% fig.PaperOrientation = 'landscape';
% saveas(gca,fullfile(Dir.figs,['DelayRespPowerTrial_glme_ACC_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']))

%%

figure('Position',[200 200 700 600]);

prc(1) = prctile(tblN_1.betaAvgD_1,25);
prc(2) = prctile(tblN_1.betaAvgD_1,75);

for gi = 1:2
    subplot(2,2,gi);hold all;axis square

    y1 = mean(tblN_1.betaAvgM(tblN_1.button_resp_corr==1 & double(tblN_1.group)==gi & tblN_1.betaAvgD_1<prc(1)));
    y0 = mean(tblN_1.betaAvgM(tblN_1.button_resp_corr==0 & double(tblN_1.group)==gi & tblN_1.betaAvgD_1<prc(1)));
    plot([0 1], [y0 y1])

    y1 = mean(tblN_1.betaAvgM(tblN_1.button_resp_corr==1 & double(tblN_1.group)==gi & tblN_1.betaAvgD_1>prc(2)));
    y0 = mean(tblN_1.betaAvgM(tblN_1.button_resp_corr==0 & double(tblN_1.group)==gi & tblN_1.betaAvgD_1>prc(2)));
    plot([0 1], [y0 y1])
    set(gca,'XLim',[-1 2],'XTick',[0 1],'XTickLabel',{'Incorrect','Correct'})
    ylabel('maintenance beta power')
    xlabel('Memory accuracy of trial N')
    title(groupStr{gi})
    legend({'low (25%) Deletion (n-1)','high (75%) Deletion (n-1)'},'Location','southoutside')
end

prc(1) = prctile(tblN_1.betaAvgM,25);
prc(2) = prctile(tblN_1.betaAvgM,75);

for gi = 1:2
    subplot(2,2,gi+2);hold all;;axis square

    y1 = mean(tblN_1.betaAvgD_1(tblN_1.button_resp_corr==1 & double(tblN_1.group)==gi & tblN_1.betaAvgM<prc(1)));
    y0 = mean(tblN_1.betaAvgD_1(tblN_1.button_resp_corr==0 & double(tblN_1.group)==gi & tblN_1.betaAvgM<prc(1)));
    plot([0 1], [y0 y1])

    y1 = mean(tblN_1.betaAvgD_1(tblN_1.button_resp_corr==1 & double(tblN_1.group)==gi & tblN_1.betaAvgM>prc(2)));
    y0 = mean(tblN_1.betaAvgD_1(tblN_1.button_resp_corr==0 & double(tblN_1.group)==gi & tblN_1.betaAvgM>prc(2)));
    plot([0 1], [y0 y1])
    set(gca,'XLim',[-1 2],'XTick',[0 1],'XTickLabel',{'Incorrect','Correct'})
    ylabel('Deletion beta power of trial n-1')
    xlabel('Memory accuracy of trial N')
    title(groupStr{gi})
    legend({'low (25%) maintenance (n)','high (75%) maintenance (n)'},'Location','southoutside')
end
% saveas(gca,fullfile(Dir.figs,['DelayRespPowerTrial_glme_ACC_MDprc',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

%%
for gi = 1:2
    tmp_tbl = tblN_1(double(tblN_1.group)==gi,:);

    %     glme = fitglme(tmp_tbl,'button_resp_corr ~ betaAvgD_1*betaAvgM*ss_num+(1+ss_num+betaAvgD_1+betaAvgM|subj)','Distribution','Binomial','Link','logit');
    glme = fitglme(tmp_tbl,'button_resp_corr ~ betaAvgD_1*ss_num+(1+ss_num+betaAvgD_1|subj)','Distribution','Binomial','Link','logit');

%     glme = fitglme(tmp_tbl,'button_resp_corr ~ betaAvgD_1*ss_num+(1|subj)','Distribution','Binomial','Link','logit');

%     t = dataset2table(glme.anova)
%     t.formula{1} = char(glme.Formula);
%     writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%         'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
%     t = dataset2table(glme.Coefficients);
%     t.formula{1} = char(glme.Formula);
%     writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%         'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
end

