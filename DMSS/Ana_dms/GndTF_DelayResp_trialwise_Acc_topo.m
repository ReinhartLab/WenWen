clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 0;% enforced for accuracy
IsLap = 1;
IsdePhase=1;
IsBL2preDelay = 1;% enforced

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Older','Both'};
condStr = {'ss1','ss2','ss4'};
ssN = [1 2 4];

freq.betaFreq = [15 25];% Hz
% freq.betaFreq = [8 12];% Hz
% freq.betaFreq = [1 3];% Hz
% freq.betaFreq = [4 7];% Hz

% freq.betaFreq = [28 40];% Hz

timeROI.Delay = [0 3];% in s
timeROI.Post = [0 0.5];% in s

%%
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

        timeID.Delay = dsearchn(tfDat.PowN.time',timeROI.Delay');
        timeID.Delay = timeID.Delay(1):timeID.Delay(2);

        M.subj = ones(height(M),1)*sub_i;
        M.group = ones(height(M),1)*subs.group(sub_i);

        for c = 1:numel(tfDat.PowN.label)
            tmp_tbl = M;

            tmp_tbl.betaAvgM = squeeze(mean(mean(tfDat.PowN.powspctrm(:,c,freqID.beta,timeID.Delay),3,"omitnan"),4,"omitnan"));
            tmp_tbl.betaAvgD = squeeze(mean(mean(tfDat.PowN_1.powspctrm(:,c,freqID.beta,timeID.Post),3,"omitnan"),4,"omitnan"));
            tmp_tbl.chan = ones(height(tmp_tbl),1)*c;
            tmp_tbl.chanName = repmat(tfDat.PowN.label(c),height(tmp_tbl),1);

            %----------------
            tmp_cidx = ismember(tmp_tbl.trl-1,tmp_tbl.trl);
            Mc = tmp_tbl(tmp_cidx,:);

            tmp_pidx = ismember(tmp_tbl.trl,Mc.trl-1);
            Mp = tmp_tbl(tmp_pidx,{'betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'});
            Mp = renamevars(Mp,{'betaAvgD','trl','button_resp_rt','button_resp_corr','ss_num'},{'betaAvgD_1','trl_1','rt_1','corr_1','ss_1'});

            Mc = [Mc Mp];
            tblN_1 = [tblN_1;Mc];
        end
    end
end

save(fullfile(Dir.results,['DelayRespBetaPow_trialwise_acc_chans',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'tblN_1','-v7.3')

%% glme for frontal beta power 
load(fullfile(Dir.results,['DelayRespBetaPow_trialwise_acc_chans',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))

%--------selected channels
tblN_1.subj = categorical(tblN_1.subj);
tblN_1.ss_num = categorical(tblN_1.ss_num); %n ss
tblN_1.ss_1 = categorical(tblN_1.ss_1);% n-1 ss

chanN = max(tblN_1.chan);
chans = unique(tblN_1.chanName,'stable');
% -------separate young and old

for gi = 1:2
    for c = 1:chanN
        c_tbl = tblN_1(tblN_1.chan==c & tblN_1.group== gi,:);
        glme = fitglme(c_tbl,'button_resp_corr ~ betaAvgD_1*ss_num+(1+betaAvgD+ss_num|subj)','Distribution','Binomial','Link','logit');

        t = dataset2table(glme.anova);
        i = 3;% main effect for deletion
        myReg.group{gi}.chans.pEEG(c) = t.pValue(i);
        myReg.group{gi}.chans.FEEG(c) = t.FStat(i);

        i = 4; % ss*deletion
        myReg.group{gi}.chans.Finter(c) = t.FStat(i);
        myReg.group{gi}.chans.pinter(c) = t.pValue(i);

        t = dataset2table(glme.Coefficients);
        i=4;
        myReg.group{gi}.chans.eEEG(c) = t.Estimate(i);
    end
end

 %% use ft to plot, so that pdf can be editable
    ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));

    cfg= [];
    cfg.latency = 0;
    cfg.channel = chans;
    ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{1:2,1});

      %%
    gi = 2;
    ft_temp.avg = myReg.group{gi}.chans.eEEG;
    figure('position',[100 100 350 350])
    cfg = [];
    cfg.colormap = jet;
    cfg.highlight          =  'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize   = 20;
    cfg.highlightchannel = chans(myReg.group{gi}.chans.pEEG<threshP);
    cfg.style = 'straight';
    cfg.marker = 'off';
    cfg.layout = 'easycapM11.mat';
    cfg.comment = 'no';
    cfg.figure = 'gca';
    ft_topoplotER(cfg,ft_temp);
 
    caxis([-0.1, 0.1])
    h = colorbar;
    h.Label.String = sprintf('GLME Coefficient (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1.2 0 0];
    h.Ticks = [-0.1:0.05:0.1];
    h.TickLength = 0;
    h.FontSize = 13;
    title(groupStr{gi},'N+1 deletion main effect');

    saveas(gcf,fullfile(Dir.figs,['DelayRespBetaPowerTopoGLME_N+1Coef_acc',groupStr{gi},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
    saveas(gcf,fullfile(Dir.figs,['DelayRespBetaPowerTopoGLME_N+1Coef_acc',groupStr{gi},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))

    %% plot B of beta

load chanLocs_dms.mat
[is, loc] = ismember(unique(tblN_1.chanName,'stable'),{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.05;
mkSize = 5;
myColors = jet;
%     myColors = hot;

figure('Position',[100 20 1000 500]);

for gi = 1:2
    subplot(2,2,gi)
    topoplot(myReg.group{gi}.chans.eEEG,chanloc,'emarker2',{find(myReg.group{gi}.chans.pEEG<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    caxis([-0.1, 0.1])
    h = colorbar;
    h.Label.String = sprintf('GLME Coef(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(groupStr{gi},'N+1 deletion main effect');

    subplot(2,2,gi+2)
    topoplot(myReg.group{gi}.chans.Finter,chanloc,'emarker2',{find(myReg.group{gi}.chans.pinter<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    caxis([0, 5])
    h = colorbar;
    h.Label.String = sprintf('GLME Coef Fvalue(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(groupStr{gi},'N+1 ss*deletion');
end

saveas(gcf,fullfile(Dir.figs,['DelayRespBetaPowerTopoGLME_N+1acc_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))