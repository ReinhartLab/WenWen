clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 0;% use acc as behavioral index
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
            dat{sub_i}.beha{cond_i} = TrlBeha{cond_i}.button_resp_corr;

            tmp_tbl = table(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.beha{cond_i},'VariableNames',{'eeg','beha'});
            tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);
            tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
            tmp_tbl.ss = ones(height(tmp_tbl),1)*cond_i;

            tbl = [tbl;tmp_tbl];
            
            for c = 1:length(tfDat{1}.label)
                tmp_pow = squeeze(dat{sub_i}.betaChans{cond_i}(:,c));
                tmp_beha =  dat{sub_i}.beha{cond_i};

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

save(fullfile(Dir.results,['DelayBetaPower_ACC',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'dat','tbl','tblChans','-v7.3')
%% glme for frontal beta power and each channel
load(fullfile(Dir.results,['DelayBetaPower_ACC',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))
tbl.subj = categorical(tbl.subj);
tbl.ss = categorical(tbl.ss);
tbl.group = categorical(tbl.group);

% first put everything as fixed effect, ss has no main effect or
% interaction with other terms (frontal chans)

glme = fitglme(tbl,'beha~eeg*group*ss+(1|subj)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group*ss+(1+ss|subj)','Distribution','Binomial','Link','logit');
glme = fitglme(tbl,'beha~eeg*group+(1|subj)+(1|ss)','Distribution','Binomial','Link','logit');
t = dataset2table(glme.anova)
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNacc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova1'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNacc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients1'])
%%
% control for ss, interaction group*eeg
glme = fitglme(tbl,'beha ~ eeg*group+(1|subj) + (1|ss)','Distribution','Binomial','Link','logit');
t = dataset2table(glme.anova);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNacc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_anova2'])
t = dataset2table(glme.Coefficients);
t.fomula{1} = char(glme.Formula);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNacc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[frontalROI{:},'_Coefficients2'])

for gi = 1:2
    tmp_tbl = tbl(double(tbl.group)==gi,:);

       glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj) + (1|ss)','Distribution','Binomial','Link','logit');
%  glme = fitglme(tmp_tbl,'beha ~ eeg*ss+(1|subj)','Distribution','Binomial','Link','logit');
    t = dataset2table(glme.anova);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNacc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
    t = dataset2table(glme.Coefficients)
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglmeNNacc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
end

for s = 1:3
    tmp_tbl = tbl(double(tbl.ss)==s,:);
    glme = fitglme(tmp_tbl,'beha ~ eeg*group+(1|subj)','Distribution','Binomial','Link','logit');
    t = dataset2table(glme.anova)
    t = dataset2table(glme.Coefficients);
end

%%



%%
tblChans.subj = categorical(tblChans.subj);
tblChans.ss = categorical(tblChans.ss);
tblChans.group = categorical(tblChans.group);

for c = 1:length(tfDat{1}.label)
    tmp_tbl = tblChans(tblChans.chans ==c,:);

    %     glme = fitglme(tmp_tbl,'beha ~ eeg*ss+(1|subj)','Distribution','Binomial','Link','logit');
    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|ss)+(1|subj)','Distribution','Binomial','Link','logit');
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

%% correlation topography beta

myFigBasic;
load chanLocs_dms.mat
[is, loc] = ismember(tfDat{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
mkSize = 8;
myColors = jet;
%% use ft to plot, so that pdf can be editable
ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));

cfg= [];
cfg.latency = 0;
cfg.channel = tfDat{1}.label;
% ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{1:2,1});

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
saveas(gcf,fullfile(Dir.figs,['DelayBetaPowerTopoGLME_NNacc_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Delay(1)),'~',num2str(timeROI.Delay(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
