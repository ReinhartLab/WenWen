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
IsBL2preDelay = 1;% baselined to pre response
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
condStr = {'ss1','ss2','ss4'};
chanROI = {'Fz','F1','F2','FCz','AFz'};
% chanROI = {'Fz','F1','F2','FCz','AFz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz

timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s
timeROI.Pre = [-0.4 0];% in s

%%
for sub_i = 15%1:subN
    subname = subs.name{sub_i};
    matName = fullfile(Dir.results,[subname,'_resp_ft_var',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);
    end
    %% time-frequency map
    myFigBasic;

    figure('Position',[100 100 900 800]);
    for cond_i = 1:3
        subplot(3,3,cond_i);
        cfg = [];
        cfg.xlim  = timeROI.all;
        cfg.ylim = [1 40];
        %         cfg.zlim = [-1 1];
        % cfg.colormap = bkr;
        cfg.colormap = jet;
        cfg.channel = chanROI;
        cfg.figure = 'gca';

        ft_singleplotTFR(cfg, tfDat{cond_i});

        title([subname,condStr{cond_i}],[chanROI{:}],'FontSize',10);
        hc = colorbar;
        hc.Label.String = 'Relative variance(a.u.)';
        hc.Label.Rotation = -90;
        hc.Label.Position = [4.2 0 0];

        set(gca,'ytick',0:10:length(tfDat{1}.freq),'XTick',timeROI.all)
        if cond_i == 1
            ylabel('\bfFrequency(Hz)');
        end
        xlabel('Time(0s=response)')
        rectangle('Position',[0 0.5 0.5 40],'LineWidth',1.2,'LineStyle','--')

        subplot(3,3,3+cond_i);
        cfg = [];
        cfg.ylim = freq.betaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        %         cfg.zlim = [-1 1];
        cfg.highlightchannel = chanROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.title = sprintf('Beta(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2));
        ft_topoplotTFR(cfg, tfDat{cond_i});colorbar

        subplot(3,3,6+cond_i);
        cfg = [];
        cfg.ylim = freq.alphaFreq;%freq
        cfg.xlim = timeROI.Post;%time
        %         cfg.zlim = [-1 1];
        cfg.highlightchannel = chanROI;
        cfg.layout = 'easyCapM1';
        cfg.colormap = jet;
        cfg.markersymbol = '.';
        cfg.highlight = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize = 18;
        cfg.figure      = 'gca';
        cfg.title = sprintf('Alpha(%d~%dHz)',freq.alphaFreq(1),freq.alphaFreq(2));

        ft_topoplotTFR(cfg, tfDat{cond_i});colorbar
    end

    saveas(gcf,fullfile(Dir.figs,'IndividualResults',[subname,'NeuVar_Resp_',[chanROI{:}],txtCell{IsLap+1,1},num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

    %%
    freqID.beta = dsearchn(tfDat{1}.freq',freq.betaFreq');
    freqID.beta = freqID.beta(1):freqID.beta(2);

    freqID.alpha = dsearchn(tfDat{1}.freq',freq.alphaFreq');
    freqID.alpha = freqID.alpha(1):freqID.alpha(2);

    timeID.all = dsearchn(tfDat{1}.time',timeROI.all');
    timeID.all = timeID.all(1):timeID.all(2);

    timeID.Post = dsearchn(tfDat{1}.time',timeROI.Post');
    timeID.Post = timeID.Post(1):timeID.Post(2);

    timeID.Pre = dsearchn(tfDat{1}.time',timeROI.Pre');
    timeID.Pre = timeID.Pre(1):timeID.Pre(2);

    clear chanID
    [log1,chanID.frontal] = ismember(chanROI,tfDat{1}.label);

    clear dat

    for cond_i = 1:3
        dat.betaAvgFrontal(:,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),1,"omitnan"));
        dat.betaAvgFrontalPre(:,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.Pre),2,"omitnan"),3,"omitnan"),1,"omitnan"));

        dat.betaCurveFrontal(:,cond_i,:) = squeeze(mean(mean(tfDat{cond_i}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.all),2,"omitnan"),1,"omitnan"));

        dat.alphaAvgFrontal(:,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i}.powspctrm(chanID.frontal(log1),freqID.alpha,timeID.Post),2,"omitnan"),3,"omitnan"),1,"omitnan"));
        dat.alphaAvgFrontalPre(:,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i}.powspctrm(chanID.frontal(log1),freqID.alpha,timeID.Pre),2,"omitnan"),3,"omitnan"),1,"omitnan"));

        dat.alphaCurveFrontal(:,cond_i,:) = squeeze(mean(mean(tfDat{cond_i}.powspctrm(chanID.frontal(log1),freqID.alpha,timeID.all),2,"omitnan"),1,"omitnan"));
    end

    %% bar & time series -beta
    addpath(genpath('D:\Toolbox\crameri_v1.08'))
    myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
    fig =figure('Position',[100 100 800 350]);

    subplot(2,2,1);hold all;axis square
    mn = dat.betaAvgFrontal;
    hb = bar(mn,'FaceAlpha',0.8);
    set(gca,'xtick',[1 2 3],'XTickLabel',condStr)
    ytickformat('%.1f')
    ylabel('Variance(a.u.)')
    title(sprintf('%s beta\nPost(%.1f~%.1fs)',[chanROI{:}],timeROI.Post(1),timeROI.Post(2)))

    subplot(2,2,3);hold all;axis square
    mn = dat.betaAvgFrontalPre;
    hb = bar(mn,'FaceAlpha',0.8);

    set(gca,'xtick',[1 2 3],'XTickLabel',condStr)
    ytickformat('%.1f')
    ylabel('Variance(a.u.)')
    title(sprintf('%s beta\nPre(%.1f~%.1fs)',[chanROI{:}],timeROI.Pre(1),timeROI.Pre(2)))

    myColors = crameri('batlowW',5);

    times = tfDat{1}.time(timeID.all);
    subplot(1,2,2);hold all;axis square
    for cond_i = 1:3
        mn =dat.betaCurveFrontal(:,cond_i);
        plot(times,mn,'color',myColors(cond_i,:))
    end

    legend(condStr,'Location','bestoutside')
    ytickformat('%.1f')
    xtickangle(0)
    ylabel('Variance(a.u.)');
    xlabel('Time(0s=response)')
    title(subname,sprintf('%s %d~%dHz',[chanROI{:}],freq.betaFreq(1),freq.betaFreq(2)))
    set(gca,'XLim',timeROI.all,'XTick',[timeROI.all(1) 0 timeROI.all(2)])
    plot(get(gca,'XLim'),[0 0],'k','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k--','HandleVisibility','off')

    saveas(gca,fullfile(Dir.figs,'IndividualResults',[subname,'RespBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
end