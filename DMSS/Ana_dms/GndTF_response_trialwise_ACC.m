clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 0; % all trials for acc
IsBL2preDelay = 1;% baselined to pre response
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Both'};
condStr = {'ss1','ss2','ss4'};

frontalROI = {'CPz','CP1','CP2','Pz','Cz'};
% frontalROI = {'FCz','Fz','F1','F2','AFz'};

freq.betaFreq = [15 25];% Hz
timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s

%%
clear subsAll
tbl = table;
tblChans = table;

for sub_i = 1:subN
    subname = subs.name{sub_i};
    disp(subname)
    matName =fullfile(Dir.results,[subname,'_resp_PowerTrialwise',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);

        freqID.beta = dsearchn(tfDat{1}.freq',freq.betaFreq');
        freqID.beta = freqID.beta(1):freqID.beta(2);

        timeID.Post = dsearchn(tfDat{1}.time',timeROI.Post');
        timeID.Post = timeID.Post(1):timeID.Post(2);

        clear chanID
        [log1,chanID.frontal] = ismember(frontalROI,tfDat{1}.label);

        for cond_i = 1:3
            dat{sub_i}.betaAvgFrontal{cond_i} = squeeze(mean(mean(mean(tfDat{cond_i}.powspctrm(:,chanID.frontal(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
            dat{sub_i}.betaChans{cond_i} = squeeze(mean(mean(tfDat{cond_i}.powspctrm(:,:,freqID.beta,timeID.Post),3,"omitnan"),4,"omitnan"));
            dat{sub_i}.beha{cond_i} = PostTrlBeha{cond_i}.button_resp_corr;
            dat{sub_i}.ss1{cond_i} = PostTrlBeha{cond_i}.ss_num;

            tmp_tbl = table(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.beha{cond_i},dat{sub_i}.ss1{cond_i},'VariableNames',{'eeg','beha','ss1'});
            tmp_tbl.ss = ones(height(tmp_tbl),1)*cond_i;
            tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
            tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);

            tbl = [tbl;tmp_tbl];

%             subsAll.behaRT(sub_i,cond_i) = mean(PostTrlBeha{cond_i}.button_resp_corr);
            for c = 1:length(tfDat{1}.label)
                tmp_pow = squeeze(dat{sub_i}.betaChans{cond_i}(:,c));
                tmp_beha =  dat{sub_i}.beha{cond_i};
                subsAll.betaChansCorr(sub_i,c,cond_i) = corr(tmp_pow,tmp_beha,'rows','pairwise');

                tmp_tbl = table(tmp_pow,tmp_beha,'VariableNames',{'eeg','beha'});
                tmp_tbl.ss = ones(height(tmp_tbl),1)*cond_i;
                tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
                tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);
                tmp_tbl.chans = ones(height(tmp_tbl),1)*c;

                tblChans = [tblChans;tmp_tbl];
            end
        end
    end
end

save(fullfile(Dir.results,['RespBetaPow_trialwiseACC',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll','dat','tbl','tblChans','-v7.3')

%% glme for frontal beta power and each channel
load(fullfile(Dir.results,['RespBetaPow_trialwiseACC',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))

%--------selected channels
tbl.subj = categorical(tbl.subj);
tbl.ss = categorical(tbl.ss);
tbl.group = categorical(tbl.group);

glme = fitglme(tbl,'beha~eeg*group*ss1+(1|subj)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group*ss1+(1+eeg+ss1|subj)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group+(1|subj)+(1|ss1)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group+(1+eeg+ss1|subj)+(1|ss1)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group+(1+eeg|subj)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group+(1+eeg|subj)+(1|ss1)','Distribution','Binomial','Link','logit');
t = dataset2table(glme.anova)
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova1'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients1'])

% control for ss, interaction group*eeg
glme = fitglme(tbl,'beha ~ eeg*group+(1|subj) + (1|ss)');
t = dataset2table(glme.anova);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova2'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients2'])

for gi = 1:2
    tmp_tbl = tbl(double(tbl.group)==gi,:);

    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj) + (1|ss)');
    t = dataset2table(glme.anova);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
end

%%
%------------loop each channel
tblChans.subj = categorical(tblChans.subj);
tblChans.ss = categorical(tblChans.ss);
tblChans.group = categorical(tblChans.group);

for c = 1:length(tfDat{1}.label)
    tmp_tbl = tblChans(tblChans.chans ==c,:);

    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|ss)+(1|subj)','Distribution','Binomial','Link','logit');

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


%% correlation topography beta

load chanLocs_dms.mat
[is, loc] = ismember(tfDat{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.01;
mkSize = 8;
myColors = jet;

figure('Name','RespPower beha corr topo','Position',[100 20 800 750]);

for cond_i = 1:3
    for gi = 1:2
        tmpsubs = find(subs.group == gi);
        [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorr(tmpsubs,:,cond_i)));
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
    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorr(tmpsubs,:,cond_i)));
    subplot(3,3,(gi-1)*3+cond_i);hold all;

    topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors  );
    caxis([-5 5])
    h = colorbar;
    h.Label.String = sprintf('T-value(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title([groupStr{gi},condStr{cond_i}],['Power@N+1 RT']);

    [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorr(tmpsubs,:,cond_i)));
end

saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoACC_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

% %% use ft to plot, so that pdf can be editable
% ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));
% 
% cfg= [];
% cfg.latency = 0;
% cfg.channel = tfDat{1}.label;
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
% cfg.highlightchannel = tfDat{1}.label(myReg.glme.chans.pEEG<threshP & myReg.glme.chans.eEEG <0);
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
[is, loc] = ismember(tfDat{1}.label,{chanloc.labels});
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

saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoGLME_N+1_ACC',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))