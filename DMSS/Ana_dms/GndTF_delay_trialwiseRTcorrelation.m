clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 1;% use RT as behavioral index
IsBL2preDelay = 0;% baselined to pre trial
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Both'};
condStr = {'ss1','ss2','ss4'};

% frontalROI = {'CPz','CP1','CP2','Pz','Cz'};
frontalROI = {'FCz','Fz','F1','F2','AFz'};

freq.betaFreq = [15 25];% Hz

timeROI.Delay = [0 3];% in s

%%
clear subsAll
tbl = table;
tblChans = table;

for sub_i = 1:subN
    subname = subs.name{sub_i};
    disp(subname)
    matName =fullfile(Dir.results,[subname,'_delay_ft_pow_trialwise',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);

        freqID.beta = dsearchn(tfDat{1}.freq',freq.betaFreq');
        freqID.beta = freqID.beta(1):freqID.beta(2);

        timeID.Delay = dsearchn(tfDat{1}.time',timeROI.Delay');
        timeID.Delay = timeID.Delay(1):timeID.Delay(2);

        clear chanID
        [log1,chanID.frontal] = ismember(frontalROI,tfDat{1}.label);

        for cond_i = 1:3
            dat{sub_i}.betaAvgFrontal{cond_i} = squeeze(mean(mean(mean(tfDat{cond_i}.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.Delay),2,"omitnan"),3,"omitnan"),4,"omitnan"));
            dat{sub_i}.betaChans{cond_i} = squeeze(mean(mean(tfDat{cond_i}.powspctrm(:,:,freqID.beta,timeID.Delay),3,"omitnan"),4,"omitnan"));
            dat{sub_i}.beha{cond_i} = TrlBeha{cond_i}.button_resp_rt;

            tmp_tbl = table(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.beha{cond_i},'VariableNames',{'eeg','beha'});
            tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);
            tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
            tmp_tbl.ss = ones(height(tmp_tbl),1)*cond_i;

            tbl = [tbl;tmp_tbl];

            [subsAll.betaFrontalcorr(sub_i,cond_i),subsAll.betaFrontal_pval(sub_i,cond_i)]= corr(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.beha{cond_i},'type','Spearman','rows','complete');

            subsAll.beha(sub_i,cond_i) = mean(TrlBeha{cond_i}.button_resp_rt);

            for c = 1:length(tfDat{1}.label)
                tmp_pow = squeeze(dat{sub_i}.betaChans{cond_i}(:,c));
                tmp_beha =  dat{sub_i}.beha{cond_i};
                subsAll.betaChansCorrRT(sub_i,c,cond_i) = corr(tmp_pow,tmp_beha,'rows','pairwise');

                tmp_tbl = table(tmp_pow,tmp_beha,'VariableNames',{'eeg','beha'});
                tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);
                tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
                tmp_tbl.ss = ones(height(tmp_tbl),1)*cond_i;
                tmp_tbl.chans = ones(height(tmp_tbl),1)*c;

                tblChans = [tblChans;tmp_tbl];
            end
        end
    end
end

save(fullfile(Dir.results,['DelayBetaPowerRT_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll','dat','tbl','tblChans','-v7.3')
%% glme for frontal beta power and each channel
load(fullfile(Dir.results,['DelayBetaPowerRT_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))
tbl.subj = categorical(tbl.subj);
tbl.ss = categorical(tbl.ss);
tbl.group = categorical(tbl.group);


glme = fitglme(tbl,'beha~eeg*group*ss+(1|subj)');
t = dataset2table(glme.anova);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNrt_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova1'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNrt_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients1'])

% control for ss, interaction group*eeg
glme = fitglme(tbl,'beha ~ eeg*group+(1|subj) + (1|ss)');
t = dataset2table(glme.anova);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNrt_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova2'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNrt_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients2'])

for gi = 1:2
    tmp_tbl = tbl(double(tbl.group)==gi,:);

    %        glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj) + (1|ss)');
    glme = fitglme(tmp_tbl,'beha ~ eeg*ss+(1|subj)');
    t = dataset2table(glme.anova)
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNrt_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNrt_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
end


for s = 1:3
    tmp_tbl = tbl(double(tbl.ss)==s,:);
    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)');
    t = dataset2table(glme.anova)
    t = dataset2table(glme.Coefficients);
end
%%
tblChans.subj = categorical(tblChans.subj);
tblChans.ss = categorical(tblChans.ss);
tblChans.group = categorical(tblChans.group);

for c = 1:length(tfDat{1}.label)
    tmp_tbl = tblChans(tblChans.chans ==c,:);

    %     glme = fitglme(tmp_tbl,'beha ~ eeg*ss+(1|subj)');
    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|ss)+(1|subj)');
    t = dataset2table(glme.anova);
    i = 2;% index for power data
    myReg.glme.chans.pEEG(c) = t.pValue(i);
    myReg.glme.chans.FEEG(c) = t.FStat(i);

    t = dataset2table(glme.Coefficients);
    myReg.glme.chans.eEEG(c) = t.Estimate(i);

    for gi = 1:2
        tmp_tbl = tblChans(double(tblChans.group)==gi & tblChans.chans ==c,:);

        glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj) + (1|ss)');
        t = dataset2table(glme.anova);
        i = 2;% index for power data
        myReg.glme.chansGroup{gi}.pEEG(c) = t.pValue(i);
        myReg.glme.chansGroup{gi}.FEEG(c) = t.FStat(i);

        t = dataset2table(glme.Coefficients);
        myReg.glme.chansGroup{gi}.eEEG(c) = t.Estimate(i);
    end
end

%% FrontalBetaRT

figure('Name','FrontalBetaRT','position',[500 500 900 900]);
for gi = 1:2
    tmpsubs = find(subs.group == gi);
    [~,p,~,stats] = ttest(subsAll.betaFrontalcorr(tmpsubs,:));
    for cond_i = 1:3
        subplot(3,3,(gi-1)*3+cond_i);hold all;axis square
        for s = 1:length(tmpsubs)
            plot(dat{tmpsubs(s)}.betaAvgFrontal{cond_i},dat{tmpsubs(s)}.beha{cond_i},'.')
        end
        title(groupStr{gi},[condStr{cond_i}])
        xlabel(['BetaPower(' [frontalROI{:}] ')'])
        ylabel('Beha(RT)')
        set(gca,'ylim',[0 1.5])

        text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalcorr(tmpsubs,cond_i)),stats.tstat(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
    end
end

gi =3;% both groups
tmpsubs = 1:height(subs);
[~,p,~,stats] = ttest(subsAll.betaFrontalcorr(tmpsubs,:));
for cond_i = 1:3
    subplot(3,3,(gi-1)*3+cond_i);hold all;axis square
    for s = 1:length(tmpsubs)
        plot(dat{tmpsubs(s)}.betaAvgFrontal{cond_i},dat{tmpsubs(s)}.beha{cond_i},'.')
        lsline
    end
    title(groupStr{gi},[condStr{cond_i}])
    xlabel(['BetaPower(' [frontalROI{:}] ')'])
    ylabel('Beha(RT)')
    set(gca,'ylim',[0 1.5])

    text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalcorr(tmpsubs,cond_i)),stats.tstat(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
end
saveas(gca,fullfile(Dir.figs,['DelayBetaPowerFrontalRT',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

%% correlation topography beta

myFigBasic;
load chanLocs_dms.mat
[is, loc] = ismember(tfDat{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.05;
mkSize = 8;
myColors = jet;

figure('Name','DelayPower beha corr topo','Position',[100 20 800 750]);

for cond_i = 1:3
    for gi = 1:2
        tmpsubs = find(subs.group == gi);
        [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
        sigChanN.betaChansCorrRT(gi,cond_i) = sum((p<threshP & stats.tstat>0));% for bar plot

        subplot(3,3,(gi-1)*3+cond_i);hold all;

        topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
        caxis([-5 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i}],['@ RT']);
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
    title([groupStr{gi},condStr{cond_i}],['Power @ RT']);

    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
end

saveas(gcf,fullfile(Dir.figs,['DelayBetaPowerTopo_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
%% use ft to plot, so that pdf can be editable
ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));

cfg= [];
cfg.latency = 0;
cfg.channel = tfDat{1}.label;
ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{1:2,1});

%%
threshP = 0.001;
figure('position',[100 100 1200 350])

subplot(1,3,1);axis square
ft_temp.avg = myReg.glme.chans.FEEG;

cfg = [];
cfg.colormap = hot;
cfg.highlight          =  'on';
cfg.highlightsymbol = '.';
cfg.highlightsize   = 20;
cfg.highlightchannel = tfDat{1}.label(myReg.glme.chans.pEEG<threshP & myReg.glme.chans.eEEG >0);
cfg.style = 'both';
cfg.marker = 'off';
cfg.layout = 'easycapM11.mat';
cfg.comment = 'no';
cfg.figure = 'gca';
ft_topoplotER(cfg,ft_temp);
%
%     topoplot(myReg.groupBoth.glme.FBeta,chanloc,'emarker2',{find(myReg.groupBoth.glme.pBeta<threshP),'*','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',flipud(hot),'style','map');
caxis([0, 20])
h = colorbar;
h.Label.String = sprintf('GLME Coefficient Ftest(p<%.3f)',threshP);
h.Label.Rotation = -90;
h.Label.Position = h.Label.Position + [1.2 0 0];
h.Ticks = [0:5:20];
h.TickLength = 0;
h.FontSize = 13;
title('N-N correlation (power*RT)','FontSize',13);

for gi = 1:2
    subplot(1,3,1+gi);axis square
    ft_temp.avg = myReg.glme.chansGroup{gi}.FEEG;

    cfg = [];
    cfg.colormap = hot;
    cfg.highlight  =  'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize   = 20;
    cfg.highlightchannel = tfDat{1}.label(myReg.glme.chansGroup{gi}.pEEG<threshP & myReg.glme.chansGroup{gi}.eEEG >0);
    cfg.style = 'both';
    cfg.marker = 'off';
    cfg.layout = 'easycapM11.mat';
    cfg.comment = 'no';
    cfg.figure = 'gca';
    ft_topoplotER(cfg,ft_temp);

    caxis([0, 20])
    h = colorbar;
    h.Label.String = sprintf('GLME Coefficient Ftest(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1.2 0 0];
    h.Ticks = [0:5:20];
    h.TickLength = 0;
    h.FontSize = 13;
    title('N-N correlation (power*RT)',groupStr{gi},'FontSize',13);
end
saveas(gcf,fullfile(Dir.figs,['DelayBetaPowerTopoGLME_NN_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
