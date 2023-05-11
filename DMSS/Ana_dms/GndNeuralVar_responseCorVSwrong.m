clear;
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 0;% must be 0
IsBL2preDelay = 1;% baselined to pre response
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Young-Old'};
condStr = {'Correct','Wrong','Correct-Wrong'};

chanROI = {'Fz','F1','F2','FCz','AFz'};% frontal cluster
% from pre-post of older group
% chanROI = {'CPz','CP1','CP2','Pz','Cz'};%centroparietal cluster

freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 12];

timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s
timeROI.Pre = [-0.4 0];% in s

iterN = 200;
%%
subsAll = cell(subN,2);
for sub_i = 1:subN
    subname = subs.name{sub_i};
    matName = fullfile(Dir.results,[subname,'_resp_var_ft_varCorrWrong',num2str(iterN),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);
        subsAll(sub_i,:) = tfDat(3,:)';
    end
end

empty_idx = cellfun(@isempty,subsAll(:,1));
subsAll(empty_idx,:) = [];
subs(empty_idx,:) = [];

for gi = 1:2
    snList = subs.group == gi;
    for c = 1:2
        gndTF{gi,c} = ft_freqgrandaverage([],subsAll{snList,c});

        cfg = [];
        cfg.keepindividual = 'yes';
        tmp = ft_freqgrandaverage(cfg,subsAll{snList,c});
        gndTF{gi,c}.indv = tmp.powspctrm;
    end
end

for gi = 1:2
    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'powspctrm';
    gndTF_SSdiff{gi} = ft_math(cfg, gndTF{gi,1},gndTF{gi,2}); %difference
end

for c = 1:2
    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'powspctrm';
    gndTF_Groupdiff{c} = ft_math(cfg, gndTF{1,c},gndTF{2,c}); %group difference
end

%% time-frequency map
myFigBasic;

w = [1 1 1];
b = [0,0,1];
k = [0,0,0];
r = [1,0,0];
bkr = createcolormap(w,b,k,r,w); % 256x3 array

figure('Position',[100 100 900 1200]);

for gi = 1:2
    for c = 1:2
        subplot(6,2,c+(gi-1)*6);
        cfg = [];
        cfg.xlim  = timeROI.all;
        cfg.ylim = [1 40];
        cfg.zlim = [-1 1];
        % cfg.colormap = bkr;
        cfg.colormap = jet;

        cfg.channel = chanROI;
        cfg.figure = 'gca';

        ft_singleplotTFR(cfg, gndTF{gi,c});

        title(condStr{c},groupStr{gi},'FontSize',14);
        hc = colorbar;
        hc.Label.String = 'Relative variance(a.u.)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(gndTF{1}.freq),'XTick',timeROI.all)
        if c == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time(0s=response)')
        rectangle('Position',[0 0.5 0.5 40],'LineWidth',1.2,'LineStyle','--')

        subplot(6,2,2+c+(gi-1)*6);
        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        cfg.zlim = [-1 1];
        cfg.highlightchannel = chanROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.title = sprintf('Beta(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2));
        ft_topoplotTFR(cfg, gndTF{gi,c});colorbar

        subplot(6,2,4+c+(gi-1)*6);
        cfg = [];
        cfg.ylim = freq.alphaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        cfg.zlim = [-1 1];
        cfg.highlightchannel = chanROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.title = sprintf('Alpha(%d~%dHz)',freq.alphaFreq(1),freq.alphaFreq(2));

        ft_topoplotTFR(cfg, gndTF{gi,c});colorbar
    end
end
saveas(gcf,fullfile(Dir.figs,['NeuVar_Resp_CorVSwrong',[chanROI{:}],txtCell{IsLap+1,1},num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%%  time frequency of group average
figure('Position',[100 100 600 600]);

cfg = [];
cfg.ylim = freq.betaFreq;%freq
cfg.xlim = timeROI.Post;%time
cfg.zlim = [-0.4 0.4];
cfg.highlightchannel = chanROI;
cfg.layout = 'easyCapM1';
cfg.colormap = jet;
cfg.markersymbol = '.';
cfg.highlight = 'on';
cfg.highlightsymbol = '.';
cfg.highlightsize = 18;
cfg.figure      = 'gca';

for c = 1:2
    subplot(2,2,c);axis square
    ft_topoplotTFR(cfg, gndTF_Groupdiff{c});
    title(condStr{c},groupStr{3});colorbar
end

for gi = 1:2
    subplot(2,2,gi+2);axis square
    ft_topoplotTFR(cfg, gndTF_SSdiff{gi});
    title(condStr{3},groupStr{gi});colorbar
end
saveas(gcf,fullfile(Dir.figs,['NeuVar_Resp_GroupDiff_CorrVSWrong',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

%%
freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID.all = dsearchn(gndTF{1}.time',timeROI.all');
timeID.all = timeID.all(1):timeID.all(2);

timeID.Post = dsearchn(gndTF{1}.time',timeROI.Post');
timeID.Post = timeID.Post(1):timeID.Post(2);

timeID.Pre = dsearchn(gndTF{1}.time',timeROI.Pre');
timeID.Pre = timeID.Pre(1):timeID.Pre(2);

clear chanID
[log1,chanID.frontal] = ismember(chanROI,gndTF{1}.label);

clear dat

for gi = 1:2
    for c = 1:2
        dat.betaAvgFrontal{gi}(:,c) = squeeze(mean(mean(mean(gndTF{gi,c}.indv(:,chanID.frontal(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.betaAvgFrontalPre{gi}(:,c) = squeeze(mean(mean(mean(gndTF{gi,c}.indv(:,chanID.frontal(log1),freqID.beta,timeID.Pre),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        dat.betaCurveFrontal{gi}(:,c,:) = squeeze(mean(mean(gndTF{gi,c}.indv(:,chanID.frontal(log1),freqID.beta,timeID.all),2,"omitnan"),3,"omitnan"));

        dat.alphaAvgFrontal{gi}(:,c) = squeeze(mean(mean(mean(gndTF{gi,c}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"));
        dat.alphaAvgFrontalPre{gi}(:,c) = squeeze(mean(mean(mean(gndTF{gi,c}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.Pre),2,"omitnan"),3,"omitnan"),4,"omitnan"));

        dat.alphaCurveFrontal{gi}(:,c,:) = squeeze(mean(mean(gndTF{gi,c}.indv(:,chanID.frontal(log1),freqID.alpha,timeID.all),2,"omitnan"),3,"omitnan"));
    end
end
%% mixed ANOVA: group*cond = 2*2
clear Xa
Xa(:,1) = [reshape(dat.betaAvgFrontal{1},size(dat.betaAvgFrontal{1},1)*size(dat.betaAvgFrontal{1},2),1);...
    reshape(dat.betaAvgFrontal{2},size(dat.betaAvgFrontal{2},1)*size(dat.betaAvgFrontal{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*2,1);ones(sum(subs.group==2)*2,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) ],...
    [1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';
[SSQs.betaFrontal, DFs.betaFrontal, MSQs.betaFrontal, Fs.betaFrontal, Ps.betaFrontal]=mixed_between_within_anova(Xa);

T = subs(:,{'name','group','groupStr'});
T.corr = [dat.betaAvgFrontal{2}(:,1);dat.betaAvgFrontal{1}(:,1)];
T.wrong = [dat.betaAvgFrontal{2}(:,2);dat.betaAvgFrontal{1}(:,2)];
writetable(T,fullfile(Dir.ana,'TableOutput',['RespBetaVar_CorrVSwrong',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']))

Xa(:,1) = [reshape(dat.betaAvgFrontalPre{1},size(dat.betaAvgFrontalPre{1},1)*size(dat.betaAvgFrontalPre{1},2),1);...
    reshape(dat.betaAvgFrontalPre{2},size(dat.betaAvgFrontalPre{2},1)*size(dat.betaAvgFrontalPre{2},2),1)];
[SSQs.betaFrontalPre, DFs.betaFrontalPre, MSQs.betaFrontalPre, Fs.betaFrontalPre, Ps.betaFrontalPre]=mixed_between_within_anova(Xa);

Xa(:,1) = [reshape(dat.alphaAvgFrontal{1},size(dat.alphaAvgFrontal{1},1)*size(dat.alphaAvgFrontal{1},2),1);...
    reshape(dat.alphaAvgFrontal{2},size(dat.alphaAvgFrontal{2},1)*size(dat.alphaAvgFrontal{2},2),1)];
[SSQs.alphaFrontal, DFs.alphaFrontal, MSQs.alphaFrontal, Fs.alphaFrontal, Ps.alphaFrontal]=mixed_between_within_anova(Xa);

Xa(:,1) = [reshape(dat.alphaAvgFrontalPre{1},size(dat.alphaAvgFrontalPre{1},1)*size(dat.alphaAvgFrontalPre{1},2),1);...
    reshape(dat.alphaAvgFrontalPre{2},size(dat.alphaAvgFrontalPre{2},1)*size(dat.alphaAvgFrontalPre{2},2),1)];
[SSQs.alphaFrontalPre, DFs.alphaFrontalPre, MSQs.alphaFrontalPre, Fs.alphaFrontalPre, Ps.alphaFrontalPre]=mixed_between_within_anova(Xa);

%% bar & time series -beta
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',2);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
fig =figure('Position',[100 100 750 350]);

subplot(1,3,3);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'FaceAlpha',0.98);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
xcord = vertcat(hb(:).XEndPoints)';
plot(xcord(1,:),dat.betaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

for gi = 1:2
    errorbar(xcord(gi,:),mn(gi,:),se(gi,:),'k','LineStyle','none','HandleVisibility','off')
end
hb = bar(mn,'FaceAlpha',0.8,'HandleVisibility','off');
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
legend(condStr,'Location','southoutside')
set(gca,'xtick',[1 2],'XTickLabel',groupStr,'XLim',[0.5 2.5])
ytickformat('%.1f')
ylabel('Variance (a.u.)')
title(sprintf('%s\nPost(%.1f~%.1fs)',[chanROI{:}],timeROI.Post(1),timeROI.Post(2)))

% plot significance
if Ps.betaFrontal{3}<.05
    for gi = 1:2
        tmpdata = dat.betaAvgFrontal{gi};
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = -0.8;
        xposi = [1 2];

        if tmp_pval<=.05
            plot(xcord(gi,xposi),[1 1]*sigH,'k','HandleVisibility','off')
            text(mean(xcord(gi,xposi)),sigH,'*','FontSize',18,'HorizontalAlignment','center')
        end
    end
end


times = gndTF{1}.time(timeID.all);
for gi = 1:2
    subplot(1,3,gi);hold all;axis square;
    for c = 1:2
        mn = squeeze(mean(dat.betaCurveFrontal{gi}(:,c,:)));
        se = squeeze(std(dat.betaCurveFrontal{gi}(:,c,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(c,:)})
    end

    legend(condStr,'Location','southoutside')
    ytickformat('%.1f')
    xtickangle(0)
    ylabel('Variability (a.u.)');
    xlabel('Time(0s=response)')
    title(groupStr{gi},sprintf('%d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)],'YLim',[-0.4 1.2])
    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k--','HandleVisibility','off')

end

saveas(gca,fullfile(Dir.figs,['RespBetaVariance_CorrVSwrong_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
saveas(gca,fullfile(Dir.figs,['RespBetaVariance_CorrVSwrong_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']))
