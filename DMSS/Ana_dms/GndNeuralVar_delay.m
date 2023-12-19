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
condStr = {'Load 1','Load 2','Load 4'};
condDiffStr = {'S2-S1','S4-S1','S4-S2'};

% chanROI = {'Cz','CP1','CP2','CPz','Pz'};
chanROI = {'Fz','F1','F2','FCz','AFz'};
% chanROI = {'POz','Oz'};

freq.betaFreq = [15 25];% Hz, beta
% freq.betaFreq = [8 12];% Hz, alpha
% freq.betaFreq = [28 40];% Hz, gamma
% freq.betaFreq = [4 7];% Hz, theta
% freq.betaFreq = [1 3];% Hz, delta,might contain nan

timeROI.all = [-5.3 4.5];% in s
SStime = {[-1.6 timeROI.all(2)],[-2.8 timeROI.all(2)],[-5.2 timeROI.all(2)]};% epoch length for different SS

timeROI.delay = [0 3];% in s
timeROI.probe = [3.0 3.6];% in s
timeROI.preDelay = [-1.2 0];% in s

%%
subsAll = cell(subN,3);
for sub_i = 1:subN
    subname = subs.name{sub_i};
    matName = fullfile(Dir.results,[subname,'_var_delay_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
%     matName = fullfile(Dir.results,[subname,'_var_delay_ft_pad20',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

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

    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'indv';
    tmp = ft_math(cfg, gndTF{gi,2},gndTF{gi,1});
    gndTF_SSdiff{gi,1}.indv = tmp.indv;

    tmp = ft_math(cfg, gndTF{gi,3},gndTF{gi,1});
    gndTF_SSdiff{gi,2}.indv = tmp.indv;

    tmp = ft_math(cfg, gndTF{gi,3},gndTF{gi,2});
    gndTF_SSdiff{gi,3}.indv = tmp.indv;

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

fig = figure('Position',[100 100 1000 800]);

for gi = 1:2
    for cond_i = 1:3
        subplot(4,3,cond_i+(gi-1)*6);

        cfg = [];
        cfg.xlim  = SStime{cond_i};
        cfg.ylim = [1 40];
        cfg.zlim = [-1 1];
        cfg.colormap = jet;
        cfg.channel = chanROI;
        cfg.figure = 'gca';

        ft_singleplotTFR(cfg, gndTF{gi,cond_i});

        title([condStr{cond_i},' ' groupStr{gi}],'FontSize',14);
        hc = colorbar;
        hc.Label.String = 'Relative variance(AU)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',[-4.8 -3.6 -2.4 -1.2 0 3])
        if cond_i == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time (0s=maintenance)')
        ylabel('Frequency(Hz)')
        rectangle('Position',[0 0.5 3 40],'LineWidth',1,'LineStyle',':')
        rectangle('Position',[-1 0.5 4 40],'LineWidth',1,'LineStyle',':')

        subplot(4,3,3+cond_i+(gi-1)*6);
        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROI.delay;%time
        cfg.zlim = [-0.6 0.6];
        cfg.highlightchannel = chanROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;

        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.comment = 'no';
        ft_topoplotTFR(cfg, gndTF{gi,cond_i});colorbar;

    end
end
saveas(gcf,fullfile(Dir.figs,['NeuVar_',[chanROI{:}],txtCell{IsLap+1,1},num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
fig.PaperOrientation = 'landscape';
saveas(fig,fullfile(Dir.figs,['NeuVar_',[chanROI{:}],txtCell{IsLap+1,1},num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))

%%  time frequency of group average

figure('Position',[100 100 800 500]);
for gi = 1:3
     subplot(2,3,3+gi);
    cfg = [];
    cfg.ylim = freq.betaFreq;%freq
    cfg.xlim = timeROI.delay;%time
        cfg.zlim = [-0.6 0.1];
        if gi ==3
              cfg.zlim = [-0.5 0.5];
  end
    cfg.highlightchannel = chanROI;
    cfg.layout = 'easyCapM1';
    cfg.colormap = jet;
    cfg.markersymbol = '.';
    cfg.highlight = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 18;
    cfg.comment = 'no';
    cfg.figure      = 'gca';
    ft_topoplotTFR(cfg, GroupAvg{gi});
   hc = colorbar;
    hc.Label.String = 'Variability(a.u.)';
    hc.Label.Rotation = -90;
    hc.Ticks = cfg.zlim;
    
    subplot(2,3,gi);axis square;hold all
cfg = [];
cfg.xlim  = [-1 5];
cfg.ylim = [1 40];
cfg.zlim = [-1 1];
% cfg.colormap = bkr;
cfg.colormap = jet;
cfg.channel = chanROI;
cfg.figure = 'gca';

    ft_singleplotTFR(cfg, GroupAvg{gi});

    title(groupStr{gi},'Average set sizes','FontSize',14);
    hc = colorbar;
    hc.Label.String = 'Variability(a.u.)';
    hc.Label.Rotation = -90;
    hc.Label.Position = [4.8 0 0];

    ylabel('\bfFrequency(Hz)')
    xlabel('\bfTime(s)')
    rectangle('Position',[0 0.5 3 40],'LineWidth',1,'LineStyle','--')
    set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',[-1 0 3],'YLim',[0.5 40])
end


saveas(gcf,fullfile(Dir.figs,['NeuVar_GroupAvg_',[cfg.channel{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
saveas(gcf,fullfile(Dir.figs,['NeuVar_GroupAvg_',[cfg.channel{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))
%% time frequency of group diff

figure('Position',[100 100 1000 400])
for cond_i = 1:3
    subplot(2,4,cond_i);axis square
    cfg = [];
    cfg.xlim  = [-1 4.5];
    cfg.ylim = [1 40];
    cfg.zlim = [-1 1];
    % cfg.colormap = bkr;
    cfg.colormap = jet;

    cfg.channel = chanROI;
    cfg.figure = 'gca';

    ft_singleplotTFR(cfg, gndTF_Groupdiff{cond_i});
    set(gca,'xtick',[-1 0 3])
    rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')
    title(condStr{cond_i},groupStr{3},'FontSize',14);colorbar

    subplot(2,4,cond_i+4);
    cfg = [];
    cfg.ylim = freq.betaFreq;%freq
    cfg.xlim = timeROI.delay;%time
    cfg.zlim = [-0.8 0.8];
    cfg.highlightchannel = chanROI;
    cfg.layout = 'easyCapM1';
    cfg.colormap = jet;
    cfg.markersymbol = '.';
    cfg.highlight = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 18;
    cfg.figure      = 'gca';
    ft_topoplotTFR(cfg, gndTF_Groupdiff{cond_i});
    ylabel('\bfFrequency(Hz)')
    xlabel('\bfTime(s)')

    title(condStr{cond_i},groupStr{3});colorbar
end

subplot(2,4,4);

cfg = [];
cfg.xlim  = [-1 4.5];
cfg.ylim = [1 40];
cfg.zlim = [-1 1];
% cfg.colormap = bkr;
cfg.colormap = jet;

cfg.channel = chanROI;
cfg.figure = 'gca';

ft_singleplotTFR(cfg, GroupAvg{3});
set(gca,'xtick',[-1 0 3])
rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')
title('Collapsed',groupStr{3},'FontSize',14);colorbar

subplot(2,4,8);
cfg = [];
cfg.ylim = freq.betaFreq;%freq
cfg.xlim = timeROI.delay;%time
cfg.zlim = [-0.8 0.8];
cfg.highlightchannel = chanROI;
cfg.layout = 'easyCapM1';
cfg.colormap = jet;
cfg.markersymbol = '.';
cfg.highlight = 'on';
cfg.highlightsymbol = '.';
cfg.highlightsize = 18;
cfg.comment = 'no';
cfg.figure      = 'gca';
ft_topoplotTFR(cfg, GroupAvg{3});
title('Maintenance',groupStr{3});colorbar

saveas(gcf,fullfile(Dir.figs,['NeuVarGroupDiff_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
saveas(gcf,fullfile(Dir.figs,['NeuVarGroupDiff_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))

% %% time frequency 2*3 interaction permutation effect
% cfg = [];
% % cfg.channel          = {'MEG', '-MLP31', '-MLO12'};
% cfg.latency          = timeROI.delay;
% % cfg.frequency        = freq.betaFreq;
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_indepsamplesF';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
% % cfg.tail             = 0;
% % cfg.clustertail      = 0;
%
%   cfg.computecritval = 'yes' ;%calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes';% calculate the p-values (default='no')
%
% cfg.alpha            = 0.05;
% cfg.numrandomization = 50;
% % prepare_neighbours determines what sensors may form clusters
% % cfg_neighb.method    = 'distance';
% % cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, dataFC);
%
% design = zeros(1,size(freqFIC_planar_cmb.powspctrm,1) + size(freqFC_planar_cmb.powspctrm,1));
% design(1,1:size(freqFIC_planar_cmb.powspctrm,1)) = 1;
% design(1,(size(freqFIC_planar_cmb.powspctrm,1)+1):(size(freqFIC_planar_cmb.powspctrm,1)+...
% size(freqFC_planar_cmb.powspctrm,1))) = 2;
%
% cfg.design           = design;
% cfg.ivar             = 1;
%
% [stat] = ft_freqstatistics(cfg, freqFIC_planar_cmb, freqFC_planar_cmb);

%% mixed ANOVA
freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

timeID.all = dsearchn(gndTF{1}.time',timeROI.all');
timeID.all = timeID.all(1):timeID.all(2);

timeID.delay = dsearchn(gndTF{1}.time',timeROI.delay');
timeID.delay = timeID.delay(1):timeID.delay(2);

timeID.probe = dsearchn(gndTF{1}.time',timeROI.probe');
timeID.probe = timeID.probe(1):timeID.probe(2);

timeID.preDelay = dsearchn(gndTF{1}.time',timeROI.preDelay');
timeID.preDelay = timeID.preDelay(1):timeID.preDelay(2);

clear chanID
[log1,chanID.frontal] = ismember(chanROI,gndTF{1}.label);

clear dat SStimeID

for gi = 1:2
    for cond_i = 1:3
        SStimeID{cond_i} = dsearchn(gndTF{1}.time',SStime{cond_i}');
        SStimeID{cond_i} = SStimeID{cond_i}(1):SStimeID{cond_i}(2);

        dat.betaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.delay),2,'omitnan'),3,'omitnan'),4,'omitnan'));
        dat.betaAvgFrontalProbe{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.probe),2),3),4));
        dat.betaAvgFrontalPreDelay{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID.preDelay),2),3),4));

        dat.betaCurveFrontal{gi}{cond_i} = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,SStimeID{cond_i}),2),3));
    end
end

% mixed ANOVA: group*cond = 2*3
clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontal{1},size(dat.betaAvgFrontal{1},1)*size(dat.betaAvgFrontal{1},2),1);...
    reshape(dat.betaAvgFrontal{2},size(dat.betaAvgFrontal{2},1)*size(dat.betaAvgFrontal{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.betaFrontal, DFs.betaFrontal, MSQs.betaFrontal, Fs.betaFrontal, Ps.betaFrontal]=mixed_between_within_anova(Xa);

T = subs(:,{'name','group','groupStr'});
T.load1 = [dat.betaAvgFrontal{2}(:,1);dat.betaAvgFrontal{1}(:,1)];
T.load2 = [dat.betaAvgFrontal{2}(:,2);dat.betaAvgFrontal{1}(:,2)];
T.load3 = [dat.betaAvgFrontal{2}(:,3);dat.betaAvgFrontal{1}(:,3)];
writetable(T,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

T.group = categorical(T.group);
Meas = table([1 2 3]','VariableNames',{'ss'});
Meas.ss = categorical(Meas.ss);
rm =  fitrm(T,'load1-load3 ~ group','WithinDesign',Meas);
ranova_tbl = ranova(rm);% main within & interaction

% restoredefaultpath; % anova(rm) sometimes won't work
% anova_tbl = anova(rm);
% 
% writetable(anova_tbl,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_ANOVA_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%         'FileType','spreadsheet','Sheet','anova');
% writetable(ranova_tbl,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_ANOVA_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
%         'FileType','spreadsheet','Sheet','ranova');

% manova(rm,'By','group')
% manova(rm)
% tbl = multcompare(rm,'Cond','ComparisonType',correctionStr);
% tbl = multcompare(rm,'group','ComparisonType',correctionStr);
% T_group = multcompare(rm,'group','By','Cond','ComparisonType',correctionStr);

Xa(:,1) = [reshape(dat.betaAvgFrontalProbe{1},size(dat.betaAvgFrontalProbe{1},1)*size(dat.betaAvgFrontalProbe{1},2),1);...
    reshape(dat.betaAvgFrontalProbe{2},size(dat.betaAvgFrontalProbe{2},1)*size(dat.betaAvgFrontalProbe{2},2),1)];
[SSQs.betaFrontalProbe, DFs.betaFrontalProbe, MSQs.betaFrontalProbe, Fs.betaFrontalProbe, Ps.betaFrontalProbe]=mixed_between_within_anova(Xa);

Xa(:,1) = [reshape(dat.betaAvgFrontalPreDelay{1},size(dat.betaAvgFrontalPreDelay{1},1)*size(dat.betaAvgFrontalPreDelay{1},2),1);...
    reshape(dat.betaAvgFrontalPreDelay{2},size(dat.betaAvgFrontalPreDelay{2},1)*size(dat.betaAvgFrontalPreDelay{2},2),1)];
[SSQs.betaFrontalPreDelay, DFs.betaFrontalPreDelay, MSQs.betaFrontalPreDelay, Fs.betaFrontalPreDelay, Ps.betaFrontalPreDelay]=mixed_between_within_anova(Xa);

%% linear fit of SS
clear line_fit line_fitProbe
for gi = 1:2
    for sub_i = 1:size(dat.betaAvgFrontal{gi},1)
        line_fit{gi}.params(sub_i,:) = polyfit([1 2 4]',dat.betaAvgFrontal{gi}(sub_i,:)',1);
%         line_fitProbe{gi}.params(sub_i,:) = polyfit([1 2 4]',dat.betaAvgFrontalProbe{gi}(sub_i,:)',1);
    end
    [~,line_fit{gi}.p,~,line_fit{gi}.Tstat] = ttest(line_fit{gi}.params(:,1));
%     [~,line_fitProbe{gi}.p,~,line_fitProbe{gi}.Tstat] = ttest(line_fitProbe{gi}.params(:,1));
end

T = subs(:,{'name','group','groupStr'});
T.slope_delay = [line_fit{2}.params(:,1);line_fit{1}.params(:,1)];
% T.slope_probe = [line_fitProbe{2}.params(:,1);line_fitProbe{1}.params(:,1)];
% writetable(T,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_Slope_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

for gi = 1:2
    fprintf('%s: mean slope = %.3f, T(%d) = %.3f, p =%.3f, d = %.3f\n',...
        groupStr{gi}, mean(line_fit{gi}.params(:,1)),line_fit{gi}.Tstat.df,line_fit{gi}.Tstat.tstat,line_fit{gi}.p,computeCohen_d(line_fit{gi}.params(:,1),line_fit{gi}.params(:,1)-line_fit{gi}.params(:,1),'paired'))
end
%% anova F value- get fcritical from finv, permutation useing cluster test for time frequency map
clear xFs
for fi = 1:numel(gndTF{2}.freq)
    for ti = 1:length(timeID.delay)
        clear tmp_dat
        for gi = 1:2
            for cond_i = 1:3
                tmp_dat{gi,cond_i} = squeeze(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal,fi,timeID.delay(ti)),2,'omitnan')); % subs*chan*freq*time
            end
        end
        clear X
        X(:,1) = [vertcat(tmp_dat{1,:});vertcat(tmp_dat{2,:})];
        X(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
        X(:,3) = [sort(repmat([1:3]',sum(subs.group==1),1));sort(repmat([1:3]',sum(subs.group==2),1))];
        X(:,4) = [repmat([1:sum(subs.group==1)]',3,1);repmat([1:sum(subs.group==2)]',3,1)+sum(subs.group==1)];

        [~, ~, ~, xFs{fi,ti}, ~]=mixed_between_within_anova(X);
    end
end

clear xFsMat
for fi = 1:numel(gndTF{2}.freq)
    for ti = 1:length(timeID.delay)
        xFsMat.inter(fi,ti) = xFs{fi,ti}{4};
        xFsMat.load(fi,ti) = xFs{fi,ti}{3};
    end
end

threshP_Finter = 0.05;

Fcritical = finv(1-threshP_Finter,2,78);% for interaction

figure('position',[200 200 800 500]);

for i = 1:2
    subplot(1,2,i);hold all;axis square
    if i ==1
        tmpdat = xFsMat.inter>Fcritical;
        title('Interaction effect binarized Fmap')
    else
        tmpdat = xFsMat.load>Fcritical;
        title('Load main effect binarized Fmap')
    end
    imagesc(gndTF{1}.time(timeID.delay),gndTF{1}.freq,tmpdat)
    % im.AlphaData = xFsMat.load>Fcritical;

    xlabel('Time(s)');ylabel('Frequency(Hz)')
    set(gca,'ylim',[1 40])

    % caxis([Fcritical, 6])
    h = colorbar;
    h.Label.String = sprintf('Ftest(p<%.3f)',threshP_Finter);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
end
saveas(gcf,fullfile(Dir.figs,['DelayBetaVarianceTimeFreqFinteraction_',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
%% bar & time series
addpath(genpath('D:\Toolbox\crameri_v1.08'))
fig = figure('Position',[100 100 1000 550]);
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);

myColors = crameri('batlowW',5);

for gi = 1:2
    subplot(3,2,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveFrontal{gi}{cond_i}));
        se = squeeze(std(dat.betaCurveFrontal{gi}{cond_i},0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(gndTF{1}.time(SStimeID{cond_i}),mn,se,{'color',myColors(cond_i,:)})
    end

    legend(condStr,'Location','northeastoutside')
    ytickformat('%.1f')
    ylabel('Variance(AU)');
    xlabel('Time (0s=maintenance)')
    title(groupStr{gi},sprintf('%s %d~%dHz',[chanROI{:}],freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XTick',[-4.8 -3.6 -2.4 -1.2 0 3 timeROI.all(2)],'YLim',[-1.1 0.5],'XLim',[-5.23 4],'YTick',-1:0.5:1)
    %     set(gca,'XTick',[-4.8 -3.6 -2.4 -1.2 0 3 4 timeROI.all(2)],'YLim',[-1.1 0.5],'XLim',[0 5],'YTick',-1:0.5:1)
    %     set(gca,'XTick',[-4.8 -3.6 -2.4 -1.2 0 3 4],'YLim',[-1.2 1.2],'XLim',[-1.4 4])

    plot(get(gca,'XLim'),[0 0],'k','LineWidth',0.5,'HandleVisibility','off')
    plot([-4.8 -4.8],get(gca,'YLim'),'Color',[1 1 1]*0.9,'LineWidth',0.5,'HandleVisibility','off')
    plot([-3.6 -3.6],get(gca,'YLim'),'Color',[1 1 1]*0.9,'LineWidth',0.5,'HandleVisibility','off')
    plot([-2.4 -2.4],get(gca,'YLim'),'Color',[1 1 1]*0.9,'LineWidth',0.5,'HandleVisibility','off')
    plot([-1.2 -1.2],get(gca,'YLim'),'Color',[1 1 1]*0.9,'LineWidth',0.5,'HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','LineWidth',0.5,'HandleVisibility','off')
    plot([3 3],get(gca,'YLim'),'k','LineWidth',0.5,'Color',[1 1 1]*0.1,'HandleVisibility','off')
    plot([3 3]+mean(beha.s1_rt(beha.group==gi)),[-2 -1.65-gi*0.1+mean(beha.s1_rt(beha.group==gi))],'color',myColors(1,:),'HandleVisibility','off')
    plot([3 3]+mean(beha.s2_rt(beha.group==gi)),[-2 -1.65-gi*0.1+mean(beha.s2_rt(beha.group==gi))],'color',myColors(2,:),'HandleVisibility','off')
    plot([3 3]+mean(beha.s3_rt(beha.group==gi)),[-2 -1.65-gi*0.1+mean(beha.s3_rt(beha.group==gi))],'color',myColors(3,:),'HandleVisibility','off')

end

subplot(2,3,4);%bar
hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

groupColor = flipud(crameri('bamako',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

for gi = 1:2
    normF = (2*gi-3);
    plot(1:3,dat.betaAvgFrontal{gi}*normF,'--','Color',[groupColor(gi,:) 0.5],'LineWidth',1,'HandleVisibility','off');
end
for gi = 1:2
    normF = (2*gi-3);
    bar(1:3,mn(gi,:)*normF,'FaceColor',groupColor(gi,:),'FaceAlpha',0.8,'BarWidth',0.4)
    if line_fit{gi}.p<.05
        errorbar(1:3,mn(gi,:)*normF,se(gi,:)*normF,'k','LineWidth',1.5,'HandleVisibility','off')
    else
        errorbar(1:3,mn(gi,:)*normF,se(gi,:)*normF,'k','LineWidth',1,'LineStyle','none','HandleVisibility','off')
    end
end

legend(groupStr,'Location','bestoutside')
set(gca,'XLim',[0.25 3.75],'xtick',[1 2 3],'XTickLabel',condStr,'YLim',[-1 1.3],'YTick',[-1.2:0.6:1.6])
yticklabels({'-1.2','-0.6','0','-0.6','-1.2'})

ylabel('Variance(AU)')
title(sprintf('%s\nbeta delay (%.1f~%.1fs)',[chanROI{:}],timeROI.delay(1),timeROI.delay(2)))

% plot significance
if Ps.betaFrontal{3}<.05
    sigH = 1;
    plot([1 3],[1 1]*sigH,'k','HandleVisibility','off')
    text(2,sigH*1.12,sprintf('%.3f',Ps.betaFrontal{3}),'FontSize',12,'HorizontalAlignment','center')
end

subplot(2,3,5);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontalProbe,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontalProbe,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'HandleVisibility','off','FaceAlpha',0);
xcord = vertcat(hb(:).XEndPoints)';

plot(xcord(1,:),dat.betaAvgFrontalProbe{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontalProbe{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);

% for gi = 1:2
%     if line_fitProbe{gi}.p<.05
%         errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineWidth',1.5,'HandleVisibility','off')
%     else
%         errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
%     end
% end

legend(condStr,'Location','bestoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')
title(sprintf('%s\nbeta probe (%.1f~%.1fs)',[chanROI{:}],timeROI.probe(1),timeROI.probe(2)))
if Ps.betaFrontalProbe{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontalProbe{gi} dat.betaAvgFrontalProbe{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-1.3,-1,3);
        xposi = {[1 2],[2 3],[1 3]};
        for ss = 1:3
            if tmp_pval(ss)<=.05
                plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
                text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
            end
        end
    end
end

subplot(2,3,6);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontalPreDelay,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontalPreDelay,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'HandleVisibility','off','FaceAlpha',0);
xcord = vertcat(hb(:).XEndPoints)';

plot(xcord(1,:),dat.betaAvgFrontalPreDelay{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontalPreDelay{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);

for gi = 1:2
    %     if line_fitPreDelay{gi}.p<.05
    %         errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineWidth',1.5,'HandleVisibility','off')
    %     else
    errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
    %     end
end
legend(condStr,'Location','bestoutside')

set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%.1f')
ylabel('Variance(AU)')
title(sprintf('%s\nbeta preDelay (%.1f~%.1fs)',[chanROI{:}],timeROI.preDelay(1),timeROI.preDelay(2)))
if Ps.betaFrontal{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontalPreDelay{gi} dat.betaAvgFrontalPreDelay{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-1.3,-1,3);
        xposi = {[1 2],[2 3],[1 3]};
        for ss = 1:3
            if tmp_pval(ss)<=.05
                plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
                text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
            end
        end
    end
end
saveas(gca,fullfile(Dir.figs,['DelayBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

fig.PaperOrientation = 'landscape';
print(fig,fullfile(Dir.figs,['DelayBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r300')
%% plot slope index

fig=figure('Position',[200 200 350 400]);
hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

groupColor = flipud(crameri('bamako',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

plot(xcord(1,:),dat.betaAvgFrontal{1},'Color',[0.8 0.8 0.8],'LineWidth',0.5,'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontal{2},'Color',[0.8 0.8 0.8],'LineWidth',0.5,'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',1);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord',mn',se','k')
legend(condStr,'Location','bestoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr,'XLim',[0.5 2.5],'YLim',[-2 0],'YTick',[-2:1:1.6],'YDir','reverse')

ylabel('Variance(AU)')
title(sprintf('%s\n%d-%dHz delay (%.1f~%.1fs)',[chanROI{:}],freq.betaFreq(1),freq.betaFreq(2),timeROI.delay(1),timeROI.delay(2)))

axes('Position',[.3 .55 .2 .15])
% axes('Position',[.75 .3 .2 .18]) % under legend
hold all;box on;
% plot(1+randn(length(line_fit{1}.params),1)*0.1,line_fit{1}.params(:,1),'.','Color',ones(3,1)*0.5)
% plot(2+randn(length(line_fit{2}.params),1)*0.1,line_fit{2}.params(:,1),'.','Color',ones(3,1)*0.5)

bar([1 2],[mean(line_fit{1}.params(:,1)) mean(line_fit{2}.params(:,1))],'w','BarWidth',0.5,'FaceAlpha',0.3)
errorbar([1 2],[mean(line_fit{1}.params(:,1)) mean(line_fit{2}.params(:,1))],[std(line_fit{1}.params(:,1))/sqrt(sum(subs.group==1)) std(line_fit{2}.params(:,1))/sqrt(sum(subs.group==2))],'k','LineStyle','none')
% set(gca,'YLim', [-0.35 0.1],'ytick',[-0.35 0],'YDir','reverse','XTick',[1 2],'XTickLabel',groupStr,'XLim',get(gca,'xlim')+[0.15 -0.15]);
set(gca,'YLim', [-0.1 0],'ytick',[-0.1 0 0.05],'YDir','reverse','XTick',[1 2],'XTickLabel',groupStr,'XLim',get(gca,'xlim')+[0.15 -0.15]);
% set(gca,'YLim', [-0.1 0.1],'ytick',[-0.1 0 0.05],'YDir','reverse','XTick',[1 2],'XTickLabel',groupStr,'XLim',get(gca,'xlim')+[0.15 -0.15]);
set(gca, 'color', 'none');
title('Slope');
for gi = 1:2
    if  line_fit{gi}.p<.05
        text(gi,-0.1,'*','HorizontalAlignment','center','FontSize',20)
    end
end
saveas(gcf,fullfile(Dir.figs,['DelayVarianceSlope_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
print(fig,fullfile(Dir.figs,['DelayVarianceSlope_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r300')

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);

behaStr = {'RT','Accuracy'};

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',5);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
myColors = [237 28 50;0 166 69;25 119 200 ]./255;
for idx = 2%1:2
    if idx == 2
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));
    else
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_rt','s2_rt','s3_rt'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_rt','s2_rt','s3_rt'}));
    end

    clear g
    for gi = 1:2
%         [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Spearman');
                    [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Pearson');

        for cond = 1:3
            Y = dat.beha{gi}(:,cond);
            X = dat.betaAvgFrontal{gi}(:,cond);

            g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
            g(gi,cond).geom_point(); %Draw a scatter plot
            g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
            g(gi,cond).set_names('x','Beta variance','y',behaStr{idx}); %Set the title of the axis
            g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
            g(gi,cond).set_title(sprintf('%s: Rho = %.3f, p = %.3f',groupStr{gi},r(cond,cond),pval(cond,cond)));
            g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
        end
        %%

        clear X Y
        Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
        X(:,1) = sort(repmat([1 2 4]',size(dat.betaAvgFrontal{gi},1),1));
        X(:,2) = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);

        mdl_124 =  fitglm(X,Y);
        mdl_124_0 = fitglm(X(:,1),Y);% without beta variability

        T = subs(subs.group ==gi,{'name','group','groupStr'});
        T = [T;T;T];% repeat three times.
        T.delayVar = X(:,2);T.ss = X(:,1);T.beha = Y;
        writetable(T,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_Regress_',groupStr{gi},behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))
        myTbl{gi} = T;
        T.group = categorical(T.group);
        T.ss = categorical(T.ss);

        clear X Y
        Y = reshape(dat.beha{gi}(:,[2 3]),size(dat.beha{gi},1)*2,1);
        X(:,1) = sort(repmat([2 4]',size(dat.betaAvgFrontal{gi},1),1));
        X(:,2) = reshape(dat.betaAvgFrontal{gi}(:,[2 3]),size(dat.betaAvgFrontal{gi},1)*2,1);
        mdl_24 =  fitglm(X,Y);
        mdl_24_0 = fitglm(X(:,1),Y);% without beta variability

        %%
        cond = 4;
        Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
        X = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);
        condArray = sort(repmat([1 2 3]',sum(subs.group==gi),1));
        g(gi,cond) = gramm('x',X,'y',Y,'color',condStr(condArray)); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,cond).geom_point(); %Draw a scatter plot
        g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
        g(gi,cond).set_names('x','Beta variance','y',behaStr{idx}); %Set the title of the axis
        %         g(gi,cond).set_color_options('chroma',3); %Set the color of the point
        g(gi,cond).set_color_options('map',myColors(1:3,:)/60,'n_color',3,'n_lightness',1,'lightness',60,'legend','merge'); %Set the color of the point
        g(gi,cond).set_title(sprintf('%s\n[1,2,4],B= %.3f, p = %.3f,Rsq= %.3f\n[2,4],B = %.3f, p = %.3f,Rsq= %.3f',groupStr{gi},mdl_124.Coefficients.Estimate(3),mdl_124.Coefficients.pValue(3),mdl_124.Rsquared.Ordinary-mdl_124_0.Rsquared.Ordinary,mdl_24.Coefficients.Estimate(3),mdl_24.Coefficients.pValue(3),mdl_24.Rsquared.Ordinary-mdl_24_0.Rsquared.Ordinary));
        g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
        if idx == 2
            g(gi,cond).axe_property('XLim',[-1.2 0.45],'YLim',[0.6 1],'YTick',[0.6:0.1:1],'XTick',[-1.2:0.4:0.8]);
        else
            g(gi,cond).axe_property('XLim',[-1.4 0.45],'YLim',[0.4 1.4],'YTick',[0.4:0.2:2],'XTick',[-1.2:0.4:0.8]);
        end

        g(gi,cond).set_color_options('map',myColors(1:3,:),'n_lightness',1);
        g(gi,cond).set_point_options('base_size',2);
        g(gi,cond).set_line_options('base_size',0.5);
    end
    fig = figure('Position',[100 100 1000 500]);
    g.draw();
    saveas(gcf,fullfile(Dir.figs,['NeuVarCorr_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
%     fig.PaperOrientation = 'landscape';
%     print(fig,fullfile(Dir.figs,['NeuVarCorr_',behaStr{idx},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']),'-dpdf','-r300','-bestfit')
end

%% fitglme to examine age*variability interaction

myTbl = [myTbl{1};myTbl{2}];
myTbl.group = categorical(myTbl.group);
myTbl.ss = categorical(myTbl.ss);
writetable(myTbl,fullfile(Dir.ana,'TableOutput',['DelayVarGLMM_YoungOld_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']))

glme = fitglme(myTbl,'beha ~ delayVar*group +(1|ss)');
% glme = fitglme(myTbl,'beha ~ delayVar*group +(1|ss)','Distribution','Gamma');
t = dataset2table(glme.anova);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_anova'])

t = dataset2table(glme.Coefficients);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_Coefficients'])

full_model = fitglme(myTbl,'beha ~ delayVar*group +(1|ss)');
reduced_model = fitglme(myTbl,'beha ~ delayVar+group +(1|ss)');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood)

full_model = fitglme(myTbl,'beha ~ delayVar*group +(1|ss)','Distribution','Gamma');
reduced_model = fitglme(myTbl,'beha ~ delayVar+group +(1|ss)','Distribution','Gamma');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood)

%% normalize accuracy distribution: accuracy was a binary variable, averaging across trials make a continous variable ranging [0 1]
% our sample's avg acc has a beta distribution; therefore, we fit a beta
% distribution to estimate a & b, then use the predicted value to run glme

myTbl.beha_Bnorm = myTbl.beha;
params = betafit(myTbl.beha);
myTbl.beha_Bnorm = betainv(myTbl.beha, params(1), params(2));
writetable(myTbl,fullfile(Dir.ana,'TableOutput',['DelayVarGLMM_YoungOld_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

% figure;histogram(myTbl.beha);hold on;histogram(myTbl.beha_Bnorm);
% xlabel('Accuracy');ylabel('Density');
% saveas(gca,fullfile(Dir.figs,'BetaNormalizationAvgACC.png'))

glme = fitglme(myTbl,'beha_Bnorm ~ delayVar*group +(1|ss)');

t = dataset2table(glme.anova)
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_anova_Bnorm'])

t = dataset2table(glme.Coefficients);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVarglme_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[[chanROI{:}],'_Coef_Bnorm'])

full_model = fitglme(myTbl,'beha_Bnorm ~ delayVar*group +(1|ss)');
reduced_model = fitglme(myTbl,'beha_Bnorm ~ delayVar+group +(1|ss)');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood)



tmp = myTbl(myTbl.group=="1",:);
full_model = fitglme(tmp,'beha_Bnorm ~ delayVar +(1|ss)');
reduced_model = fitglme(tmp,'beha_Bnorm ~ 1+(1|ss)');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood)


%% use fitglme to examine age*ss interaction on variability, should be consistent with Mixed ANOVA section
glme = fitglme(myTbl,'delayVar ~ ss*group +(1|name)');
t = dataset2table(glme.anova);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_ANOVA_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_anova_glme'])

t = dataset2table(glme.Coefficients);
writetable(t,fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_ANOVA_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
    'FileType','spreadsheet','Sheet',[chanROI{:},'_Coef_glme'])

full_model = fitglme(myTbl,'delayVar ~ ss*group +(1|name)');
reduced_model = fitglme(myTbl,'delayVar ~ ss+group +(1|name)');
bayes_factor = exp(full_model.LogLikelihood- reduced_model.LogLikelihood)

return
%% transition time point: find group average, define range = [-0.1 0.1], search individual peak time within this range
% detectRange = [3.1 3.5];
% clear avgPeak
% for gi = 1:2
%     for cond_i = 1:3
%         tmp_series = mean(dat.betaCurveFrontal{gi}{cond_i});
%         [pks,locs] = findpeaks(-tmp_series);
%         tmpID = gndTF{1}.time(SStimeID{cond_i}(locs))<=detectRange(2) & gndTF{1}.time(SStimeID{cond_i}(locs))>detectRange(1);
%
%         avgPeak.pks(gi,cond_i) = -pks(tmpID);
%         avgPeak.locs(gi,cond_i) = gndTF{1}.time(SStimeID{cond_i}(locs(tmpID)));
%
%     end
% end

%% Delay Avg vs Probe Avg

clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontal{1},size(dat.betaAvgFrontal{1},1)*size(dat.betaAvgFrontal{1},2),1);...
    reshape(dat.betaAvgFrontalProbe{1},size(dat.betaAvgFrontalProbe{1},1)*size(dat.betaAvgFrontalProbe{1},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==1)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)]]';

stats = rm_anova2(Xa(:,1),Xa(:,4),Xa(:,2),Xa(:,3),{'Time','SS'})

mean(dat.betaAvgFrontalProbe{1})- mean(dat.betaAvgFrontal{1})+1

%% correlation between indvidual SS slope and baseline index (SS1 beha or variance)
figure('position',[100 100 700 1000]);

for gi = 1:2
    SS_slope = line_fit{gi}.params(:,1);

    [SSslopeR.r.var{gi},SSslopeR.p.var{gi}] = corr(dat.betaAvgFrontal{gi}, SS_slope,'type','Spearman');
    [SSslopeR.r.beha{gi},SSslopeR.p.beha{gi}] = corr(dat.beha{gi}, SS_slope,'type','Spearman');

    for cond_i = 1:3
        subplot(4,3,(gi-1)*3+cond_i);
        scatter(dat.betaAvgFrontal{gi}(:,cond_i),SS_slope)
        lsline
        xlabel('Delay variance')
        ylabel('Delay slope')

        title([groupStr{gi},' ',condStr{cond_i},' var'],sprintf('Rho = %.3f, p = %.3f',SSslopeR.r.var{gi}(cond_i),SSslopeR.p.var{gi}(cond_i)))
    end

    for cond_i = 1:3
        subplot(4,3,(gi-1)*3+cond_i+6);
        scatter(dat.beha{gi}(:,cond_i),SS_slope,'green')
        lsline
        xlabel('Behav')
        ylabel('Delay slope')
        title([groupStr{gi},' ',condStr{cond_i},' behav'],sprintf('Rho = %.3f, p = %.3f',SSslopeR.r.beha{gi}(cond_i),SSslopeR.p.beha{gi}(cond_i)))
    end
end
saveas(gcf,fullfile(Dir.figs,['DelayVarCorrSlope_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%% delay load effect topography

myFigBasic;

load chanLocs_dms.mat
[is, loc] = ismember(gndTF{1}.label,{chanloc.labels});
chanloc  = chanloc(loc(is));
threshP = 0.05;
mkSize = 3;
% myColors = jet;
myColors = flipud(hot(60));
myColors = myColors(1:50,:);

% myColors = crameri('vik',256);
% myColors = myColors(256/2:256,:);

figure('Name','Load-sensitive chans','Position',[100 200 2400 900]);
timeROI.bins = {[-1 -0.5];[-0.5 0];[0 0.5];[0.5 1];[1 1.5];[1.5 2];[2 3];[-1 0];[1 2];[2 3];[0 3]};

for t = 1:length(timeROI.bins) % loop time roi
    tmpID = dsearchn(gndTF{1}.time',timeROI.bins{t}');
    tmpID = tmpID(1):tmpID(2);

    for gi = 1:2
        clear X
        for c = 1:length(gndTF{1}.label)

            load1 = squeeze(mean(mean(gndTF{gi,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
            load2 = squeeze(mean(mean(gndTF{gi,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
            load3 = squeeze(mean(mean(gndTF{gi,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN

            X(:,1) = [load1;load2;load3];
            X(:,2) = sort(repmat([1 2 3]',sum(subs.group==gi),1));
            X(:,3) = repmat([1:sum(subs.group==gi)]',3,1);

            [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
        end
        subplot(4,length(timeROI.bins),(gi-1)*length(timeROI.bins)+t);hold all

        topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
        caxis([0, 6])
        h = colorbar;
        h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.bins{t}(1),timeROI.bins{t}(2)),['N=' num2str(sum(subs.group==gi))]);
    end

    %     tmp_groupDiff = abs(diff(Ftest.p(1:2,:)<threshP,1,1));% old | young
    %     tmp_groupDiff = diff(Ftest.p(1:2,:)<threshP,1,1)==1;% old
    %     tmp_groupDiff = diff(Ftest.p(1:2,:)<threshP,1,1)==-1;% young

    gi = 3;
    clear X
    for c = 1:length(gndTF{1}.label)

        load1 = squeeze(mean(mean(gndTF{1,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load1 = [load1;tmp];

        load2 = squeeze(mean(mean(gndTF{1,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load2 = [load2;tmp];

        load3 = squeeze(mean(mean(gndTF{1,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load3 = [load3;tmp];

        X(:,1) = [load1;load2;load3];
        X(:,2) = sort(repmat([1 2 3]',subN,1));
        X(:,3) = repmat([1:subN]',3,1);

        [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
    end

    subplot(4,length(timeROI.bins),3*length(timeROI.bins)+t);hold all
    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
    caxis([0, 6])
    h = colorbar;
    %             set(h,'ytick',[0.33:0.03:0.4])
    h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(sprintf('All:%.1f~%.1fs', timeROI.bins{t}(1),timeROI.bins{t}(2)));

    subplot(4,length(timeROI.bins),2*length(timeROI.bins)+t);hold all
    %     topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(tmp_groupDiff),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
    %     caxis([0, 6])
    %     h = colorbar;
    %     h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
    %     h.Label.Rotation = -90;
    %     h.Label.Position = h.Label.Position + [1 0 0];
    %     title(sprintf('Diff:%.1f~%.1fs', timeROI.bins{t}(1),timeROI.bins{t}(2)));

    % mixed ANOVA: group*cond = 2*3

    clear tFs tPs
    for c = 1:length(gndTF{1}.label)

        load1 = squeeze(mean(mean(gndTF{1,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load1 = [load1;tmp];

        load2 = squeeze(mean(mean(gndTF{1,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load2 = [load2;tmp];

        load3 = squeeze(mean(mean(gndTF{1,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load3 = [load3;tmp];

        clear X
        X(:,1) = [load1;load2;load3];
        X(:,3:4) = fliplr(fullfact([subN 3]));
        X(:,2) = flipud([subs.group;subs.group;subs.group]);

        [~, ~, ~, tFs{c}, tPs{c}]=mixed_between_within_anova(X);
    end
    p_vals = cell2mat([tPs{:}]);
    F_vals = cell2mat([tFs{:}]);
    topoplot(F_vals(3,:),chanloc,'emarker2',{find(p_vals(3,:)<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
    caxis([0, 6])
    h = colorbar;
    h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(sprintf('2*3interaction\n%.1f~%.1fs', timeROI.bins{t}(1),timeROI.bins{t}(2)));

end
saveas(gcf,fullfile(Dir.figs,['DelayVarTopo_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%%

figure('Name','Load-sensitive chans','Position',[0 0 600 3000]);
timeROI.bins = {[-1 -0.5];[-0.5 0];[0 0.5];[0.5 1];[1 1.5];[1.5 2];[2 3];[-1 0];[1 2];[2 3];[0 3]};

for t = length(timeROI.bins)%1:length(timeROI.bins) % loop time roi
    tmpID = dsearchn(gndTF{1}.time',timeROI.bins{t}');
    tmpID = tmpID(1):tmpID(2);

    for gi = 1:2
        clear X
        for c = 1:length(gndTF{1}.label)

            load1 = squeeze(mean(mean(gndTF{gi,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
            load2 = squeeze(mean(mean(gndTF{gi,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
            load3 = squeeze(mean(mean(gndTF{gi,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN

            X(:,1) = [load1;load2;load3];
            X(:,2) = sort(repmat([1 2 3]',sum(subs.group==gi),1));
            X(:,3) = repmat([1:sum(subs.group==gi)]',3,1);

            [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
        end
        subplot(length(timeROI.bins),3,gi+(t-1)*3);hold all

        topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
        caxis([0, 6])
        h = colorbar;
        h.Label.String = sprintf('Ftest(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title(sprintf('%s:%.1f~%.1fs', groupStr{gi},timeROI.bins{t}(1),timeROI.bins{t}(2)),['N=' num2str(sum(subs.group==gi))]);
    end

    clear X
    for c = 1:length(gndTF{1}.label)

        load1 = squeeze(mean(mean(gndTF{1,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,1}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load1 = [load1;tmp];

        load2 = squeeze(mean(mean(gndTF{1,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,2}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load2 = [load2;tmp];

        load3 = squeeze(mean(mean(gndTF{1,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        tmp = squeeze(mean(mean(gndTF{2,3}.indv(:,c,freqID.beta,tmpID),3),4));%subN
        load3 = [load3;tmp];

        X(:,1) = [load1;load2;load3];
        X(:,2) = sort(repmat([1 2 3]',subN,1));
        X(:,3) = repmat([1:subN]',3,1);

        [Ftest.F(gi,c),Ftest.p(gi,c)] = RMAOV1(X);
    end

    subplot(length(timeROI.bins),3,3*t);hold all
    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
    caxis([0, 6])
    h = colorbar;
    %             set(h,'ytick',[0.33:0.03:0.4])
    h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(sprintf('Both:%.1f~%.1fs\nLoad-sensitive', timeROI.bins{t}(1),timeROI.bins{t}(2)));
end
saveas(gcf,fullfile(Dir.figs,['DelayVarTopoHorizontal_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%%
for t = 1:length(timeROI.bins) % loop time roi
    tmpID = dsearchn(gndTF{1}.time',timeROI.bins{t}');
    tmpID = tmpID(1):tmpID(2);

    fig = figure('Name','Load-sensitive chans','Position',[100 200 300 300]);
    hold all;
    topoplot(Ftest.F(gi,:),chanloc,'emarker2',{find(Ftest.p(gi,:)<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors);
    caxis([0, 8])
    h = colorbar;
    %             set(h,'ytick',[0.33:0.03:0.4])
    h.Label.String = sprintf('Ftest (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title(sprintf('All:%.1f~%.1fs', timeROI.bins{t}(1),timeROI.bins{t}(2)));
    print(fig,fullfile(Dir.figs,['DelayVarTopoJustAll_',num2str(t),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']),'-dpdf','-r300')
end