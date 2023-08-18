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
IsLap = 1;
IsdePhase=1;
IsBL2preDelay = 1;% enforced

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Both'};
condStr = {'ss1','ss2','ss4'};

% chanROI = {'CPz','CP1','CP2','Pz','Cz'};
chanROI = {'FCz','Fz','F1','F2','AFz'};

freq.betaFreq = [15 25];% Hz
% freq.alphaFreq = [8 12];% Hz

timeROI.Delay = [0 3];% in s
timeROI.Post = [0 0.5];% in s

%%
tbl  = table; % for maintenance analysis,table contains RT(n), delay power(n)
tblN_1 = table;% for n-1 analysis, table contains RT(n), delay power(n), deletion power(n-1)
tblP_1 = table; % for n+1 analysis, table contains RT(n+1), deletion power(n)
for sub_i = 1:subN
    subname = subs.name{sub_i};
    disp(subname)
    matName =fullfile(Dir.results,[subname,'_delayresp_trialwise',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);

        freqID.beta = dsearchn(tfDat.PowN.freq',freq.betaFreq');
        freqID.beta = freqID.beta(1):freqID.beta(2);

        timeID.Post = dsearchn(tfDat.PowN.time',timeROI.Post');
        timeID.Post = timeID.Post(1):timeID.Post(2);

        timeID.Delay = dsearchn(tfDat.PowN.time',timeROI.Delay');
        timeID.Delay = timeID.Delay(1):timeID.Delay(2);

        clear chanID
        [log1,chanID.frontal] = ismember(chanROI,tfDat.PowN.label);

        M.subj = ones(height(M),1)*sub_i;
        M.group = ones(height(M),1)*subs.group(sub_i);

        M.betaAvgM = squeeze(mean(mean(mean(tfDat.PowN.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.Delay),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        M.betaAvgD = squeeze(mean(mean(mean(tfDat.PowN_1.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        tbl = [tbl;M];%

        %----------------
        tmp_cidx = ismember(M.trl-1,M.trl);
        Mc = M(tmp_cidx,:);

        tmp_pidx = ismember(M.trl,Mc.trl-1);
        Mp = M(tmp_pidx,{'betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'});
        Mp = renamevars(Mp,{'betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'},{'betaAvgD_1','trl_1','rt_1','corr_1','ss_1'});

        Mc = [Mc Mp];
        tblN_1 = [tblN_1;Mc];

        %----------------
        tmp_cidx = ismember(M.trl+1,M.trl);
        Mc = M(tmp_cidx,{'subj','group','betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'});

        tmp_pidx = ismember(M.trl,Mc.trl+1);
        Mp = M(tmp_pidx,{'trl','button_resp_rt','button_resp_corr','ss_num'});
        Mp = renamevars(Mp,{'trl','button_resp_rt','button_resp_corr','ss_num'},{'trl_1','rt_1','corr_1','ss_1'});

        Mc = [Mc Mp];
        tblP_1 = [tblP_1;Mc];
    end
end

save(fullfile(Dir.results,['DelayRespBetaPow_trialwise',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'tbl','tblN_1','-v7.3')

%% glme for frontal beta power and each channel
load(fullfile(Dir.results,['DelayRespBetaPow_trialwise',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))

%--------selected channels
tblN_1.subj = categorical(tblN_1.subj);
tblN_1.ss_num = categorical(tblN_1.ss_num);
tblN_1.group = categorical(tblN_1.group);
tblN_1.ss_1 = categorical(tblN_1.ss_1);

% first put everything as fixed effect

% glme = fitglme(tblN_1,'button_resp_rt ~ betaAvgM*betaAvgD_1*group*ss_num+(1|subj)');
glme = fitglme(tblN_1,'button_resp_rt ~ betaAvgM*betaAvgD_1*group+(1|subj)');
% glme = fitglme(tblN_1,'button_resp_rt ~ betaAvgM*betaAvgD_1*group+(1|subj)+(1|ss_num)');
% glme = fitglme(tblN_1,'button_resp_rt ~ betaAvgM*betaAvgD_1*group+(1|subj)+(1|ss_num)+(1|ss_1)');

t = dataset2table(glme.anova)
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_anova1'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_Coefficients1'])
%%
for gi = 1:2
    tmp_tbl = tblN_1(double(tblN_1.group)==gi,:);

%     glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM*ss_num+(1|subj)');
%         glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM*ss_num*betaAvgD_1+(1|subj)');

        glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM*betaAvgD_1+(1|subj)');
%         glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM*betaAvgD_1+(1|subj)+(1|ss_num)');
%         glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM*betaAvgD_1+(1|subj)+(1|ss_num)+(1|ss_1)');
t = dataset2table(glme.anova)
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
end
%%

for gi = 1%:2
    for s = 1:3
    tmp_tbl = tblN_1(double(tblN_1.group)==gi & double(tblN_1.ss_num)==s,:);

        glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM*betaAvgD_1+(1|subj)');

t = dataset2table(glme.anova)

    t.fomula{1} = char(glme.Formula);
%     writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%         'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
%     writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%         'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
    end
end
%% maintenance
load(fullfile(Dir.results,['DelayRespBetaPow_trialwise',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))

%--------selected channels
tbl.subj = categorical(tbl.subj);
tbl.ss_num = categorical(tbl.ss_num);
tbl.group = categorical(tbl.group);

% first put everything as fixed effect

glme = fitglme(tbl,'button_resp_rt ~ betaAvgM*group+(1|subj)+(1|ss_num)');
t = dataset2table(glme.anova)
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_delay_a1'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_delay_Coef1'])
%%
for gi = 1:2
    tmp_tbl = tbl(double(tbl.group)==gi,:);

    glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM+(1|subj)+(1|ss_num)');
    t = dataset2table(glme.anova);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_delay_a1'])
    t = dataset2table(glme.Coefficients)
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_delay_Coef1'])
end

%% deletion
load(fullfile(Dir.results,['DelayRespBetaPow_trialwise',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))

%--------selected channels
tblP_1.subj = categorical(tblP_1.subj);
tblP_1.group = categorical(tblP_1.group);
tblP_1.ss_1 = categorical(tblP_1.ss_1);

% first put everything as fixed effect

glme = fitglme(tblP_1,'rt_1 ~ betaAvgD*group+(1|subj)+(1|ss_1)');
t = dataset2table(glme.anova)
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_resp_a1'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_resp_Coef1'])
%%
for gi = 1:2
    tmp_tbl = tbl(double(tbl.group)==gi,:);

    glme = fitglme(tmp_tbl,'button_resp_rt ~ betaAvgM+(1|subj)+(1|ss_num)');
    t = dataset2table(glme.anova);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_delay_a1'])
    t = dataset2table(glme.Coefficients)
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayRespPowerTrial_glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_delay_Coef1'])
end
%%
%------------loop each channel
tblChans.subj = categorical(tblChans.subj);
tblChans.ss = categorical(tblChans.ss);
tblChans.group = categorical(tblChans.group);

for c = 1:length(tfDat.PowN.label)
    tmp_tbl = tblChans(tblChans.chans ==c,:);

    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|ss)+(1|subj)');

    t = dataset2table(glme.anova);
    i = 2;% index for power data
    myReg.glme.chans.pEEG(c) = t.pValue(i);
    myReg.glme.chans.FEEG(c) = t.FStat(i);

    t = dataset2table(glme.Coefficients);
    myReg.glme.chans.eEEG(c) = t.Estimate(i);

       % -------separate young and old

    for gi = 1:2
        tmp_tbl = tblChans(tblChans.chans ==c & double(tblChans.group)== gi,:);

        glme = fitglme(tmp_tbl,'beha ~ 1+eeg +(1|subj) +(1|ss)');
        i = 2;% index for power data
        t = dataset2table(glme.anova);
        myReg.group{gi}.glme.pBeta(c) = t.pValue(i);
        myReg.group{gi}.glme.FBeta(c) = t.FStat(i);

        t = dataset2table(glme.Coefficients);
        myReg.group{gi}.glme.eBeta(c) = t.Estimate(i);
    end
end

%% FrontalBetaRT

figure('Name','FrontalBetaRT','position',[500 500 900 900]);
for gi = 1:2
    tmpsubs = find(subs.group == gi);
    [~,p,~,stats] = ttest(subsAll.betaFrontalRTcorr(tmpsubs,:));
    for cond_i = 1:3
        subplot(3,3,(gi-1)*3+cond_i);hold all;axis square
        for s = 1:length(tmpsubs)
            plot(dat{tmpsubs(s)}.betaAvgFrontal{cond_i},dat{tmpsubs(s)}.behaRT{cond_i},'.')
        end
        title(groupStr{gi},[condStr{cond_i}])
        xlabel(['BetaPower(' [chanROI{:}] ')'])
        ylabel('Beha(RT)')
        set(gca,'ylim',[0 1.5])

        text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalRTcorr(tmpsubs,cond_i)),stats.tstat(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalRTcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
    end
end

gi =3;% both groups
tmpsubs = 1:height(subs);
[~,p,~,stats] = ttest(subsAll.betaFrontalRTcorr(tmpsubs,:));
for cond_i = 1:3
    subplot(3,3,(gi-1)*3+cond_i);hold all;axis square
    for s = 1:length(tmpsubs)
        plot(dat{tmpsubs(s)}.betaAvgFrontal{cond_i},dat{tmpsubs(s)}.behaRT{cond_i},'.')
        lsline
    end
    title(groupStr{gi},[condStr{cond_i}])
    xlabel(['BetaPower(' [chanROI{:}] ')'])
    ylabel('Beha(RT)')
    set(gca,'ylim',[0 1.5])

    text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalRTcorr(tmpsubs,cond_i)),stats.tstat(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalRTcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
end
% saveas(gca,fullfile(Dir.figs,['RespBetaPowerFrontalRT',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

%% correlation topography beta

load chanLocs_dms.mat
[is, loc] = ismember(tfDat.PowN.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.05;
mkSize = 8;
myColors = jet;

figure('Name','RespPower beha corr topo','Position',[100 20 800 750]);

for cond_i = 1:3
    for gi = 1:2
        tmpsubs = find(subs.group == gi);
        [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
        sigChanN.betaChansCorrRT(gi,cond_i) = sum((p<threshP & stats.tstat<0));% for bar plot

        subplot(3,3,(gi-1)*3+cond_i);hold all;

        topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
        caxis([-5 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i}],['Power@N+1 RT']);
    end

    gi =3;% both
    tmpsubs = 1:height(subs);
    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
    subplot(3,3,(gi-1)*3+cond_i);hold all;

    topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors  );
    caxis([-5 5])
    h = colorbar;
    h.Label.String = sprintf('T-value(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title([groupStr{gi},condStr{cond_i}],['Power@N+1 RT']);

    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
end

% saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopo_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%% split young and old based on their correlatoin coefficients
myColors = jet;
cond_i = 3;
threshP = 0.01;
figure('position',[200 200 600 600]);
for gi = 1:2
    tmpsubs = find(subs.group == gi);
    tmpBeha = subsAll.behaRT(tmpsubs,cond_i);

    tmp_upper = tmpBeha>median(tmpBeha);

    subplot(2,2,(gi-1)*2+1)
    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs(tmp_upper),:,cond_i)));
    topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    caxis([-5 5])
    h = colorbar;
    h.Label.String = sprintf('T-value(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title([groupStr{gi},condStr{cond_i},'-highRT'],sprintf('Power*n+1 RT, subjN=%d',sum(tmp_upper)));

    subplot(2,2,(gi-1)*2+2)
    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs(~tmp_upper),:,cond_i)));
    topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    caxis([-5 5])
    h = colorbar;
    h.Label.String = sprintf('T-value(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title([groupStr{gi},condStr{cond_i},'-lowRT'],sprintf('Power*n+1 RT, subjN=%d',sum(~tmp_upper)));
end

% saveas(gcf,fullfile(Dir.figs,['RespBetaPowerN+1_RT_TopoSplitRT',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

% %% use ft to plot, so that pdf can be editable
% ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));
% 
% cfg= [];
% cfg.latency = 0;
% cfg.channel = tfDat.PowN.label;
% ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{1:2,1});
% 
% %%
% ft_temp.avg = myReg.glme.chans.FEEG;
% threshP = 0.05;
% figure('position',[100 100 350 350])
% cfg = [];
% cfg.colormap = hot;
% cfg.highlight          =  'on';
% cfg.highlightsymbol = '.';
% cfg.highlightsize   = 20;
% cfg.highlightchannel = tfDat.PowN.label(myReg.glme.chans.pEEG<threshP & myReg.glme.chans.eEEG <0);
% cfg.style = 'both';
% cfg.marker = 'off';
% cfg.layout = 'easycapM11.mat';
% cfg.comment = 'no';
% cfg.figure = 'gca';
% ft_topoplotER(cfg,ft_temp);
% 
% caxis([0, 6])
% h = colorbar;
% h.Label.String = sprintf('GLME Coefficient Ftest(p<%.3f)',threshP);
% h.Label.Rotation = -90;
% h.Label.Position = h.Label.Position + [1.2 0 0];
% h.Ticks = [0:2:6];
% h.TickLength = 0;
% h.FontSize = 13;
% title('N+1 correlation (power*RT)','FontSize',13);
% 
% saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoGLME_N+1_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%% plot B of beta

load chanLocs_dms.mat
[is, loc] = ismember(tfDat.PowN.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.01;
mkSize = 5;
myColors = jet;
%     myColors = hot;
%     myColors = myColors(256/2:256,:);

figure('Position',[100 20 1000 500]);

for gi = 1:2
    subplot(1,3,gi)
    topoplot(myReg.group{gi}.glme.FBeta,chanloc,'emarker2',{find(myReg.group{gi}.glme.pBeta<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    caxis([0, 5])
    h = colorbar;
    h.Label.String = sprintf('GLME Coef Fvalue(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(groupStr{gi},'N+1 corrRT');
end
subplot(1,3,3)
topoplot(myReg.glme.chans.FEEG,chanloc,'emarker2',{find(myReg.glme.chans.pEEG<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
caxis([0, 5])
h = colorbar;
h.Label.String = sprintf('GLME Coef Fvalue(p<%.3f)',threshP);
h.Label.Rotation = -90;
h.Label.Position = h.Label.Position + [1 0 0];
title('Both: eeg','N+1 corr RT');

% saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoGLME_N+1_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))