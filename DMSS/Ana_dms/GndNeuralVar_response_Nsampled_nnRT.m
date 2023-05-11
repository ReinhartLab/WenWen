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

groupStr = {'Young','Old','Both'};
condStr = {'ss1','ss2','ss4'};

% frontalROI = {'CPz','CP1','CP2','Pz','Cz'};
frontalROI = {'FCz','Fz','F1','F2','AFz'};

% occipROI = {'POz','Oz'};
freq.betaFreq = [15 25];% Hz
% freq.alphaFreq = [8 13];% Hz

timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s
timeROI.Pre = [-0.4 0];% in s

%%
for Nsamps = 12%[6:2:12]
    clear subsAll
    for sub_i = 1:subN
        subname = subs.name{sub_i};
        matName =fullfile(Dir.results,[subname,'_resp_ft_var_nnCorr_NsampledAT',num2str(Nsamps),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

        if isfile(matName)
            load(matName);

            freqID.beta = dsearchn(tfDat{1}{1}.freq',freq.betaFreq');
            freqID.beta = freqID.beta(1):freqID.beta(2);

%             freqID.alpha = dsearchn(tfDat{1}{1}.freq',freq.alphaFreq');
%             freqID.alpha = freqID.alpha(1):freqID.alpha(2);

            timeID.Post = dsearchn(tfDat{1}{1}.time',timeROI.Post');
            timeID.Post = timeID.Post(1):timeID.Post(2);

            clear chanID
            [log1,chanID.frontal] = ismember(frontalROI,tfDat{1}{1}.label);
%             [log2,chanID.occip] = ismember(occipROI,tfDat{1}{1}.label);

            for cond_i = 1:3
                for b = 1:length(tfDat{cond_i})
                    dat{sub_i}.betaAvgFrontal{cond_i}(b) = squeeze(mean(mean(mean(tfDat{cond_i}{b}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));
                    dat{sub_i}.betaChans{cond_i}(b,:) = squeeze(mean(mean(tfDat{cond_i}{b}.powspctrm(:,freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"));
                    dat{sub_i}.behaRT{cond_i}(b) = mean(TrlBeha{cond_i}{b}.button_resp_rt,'omitnan');
                end

                subsAll.betaFrontalVar(sub_i,cond_i) = mean(dat{sub_i}.betaAvgFrontal{cond_i},'omitnan');
                subsAll.behaRT(sub_i,cond_i) = mean(dat{sub_i}.behaRT{cond_i},'omitnan');

                [subsAll.betaFrontalRTcorr(sub_i,cond_i),subsAll.betaFrontalRT_pval(sub_i,cond_i)]= corr(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.behaRT{cond_i},'type','Spearman');
    
                 for c = 1:length(tfDat{1}{1}.label)
                    tmp_var = squeeze(dat{sub_i}.betaChans{cond_i}(:,c));

                    tmp_beha =  dat{sub_i}.behaRT{cond_i};
                    subsAll.betaChansCorrRT(sub_i,c,cond_i) = corr(tmp_var,tmp_beha');
                 end
            end         
        end
    end
    save(fullfile(Dir.results,['RespBetaVar_nnRTSampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll','-v7.3')

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
            title(groupStr{gi},[condStr{cond_i},' sampN=' num2str(Nsamps)])
            xlabel(['BetaVariance(' [frontalROI{:}] ')'])
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
        title(groupStr{gi},[condStr{cond_i},' sampN=' num2str(Nsamps)])
        xlabel(['BetaVariance(' [frontalROI{:}] ')'])
        ylabel('Beha(RT)')
        set(gca,'ylim',[0 1.5])

        text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalRTcorr(tmpsubs,cond_i)),stats.tstat(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalRTcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
    end
    saveas(gca,fullfile(Dir.figs,['RespBetaVarianceFrontalRT_nnSampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

   
    %% correlation topography beta

    myFigBasic;

    load chanLocs_dms.mat
    [is, loc] = ismember(tfDat{1}{1}.label,{chanloc.labels});
    chanloc  = chanloc(loc(is));
    threshP = 0.01;
    mkSize = 8;
    myColors = jet;

    % myColors = crameri('vik',256);
    % myColors = myColors(256/2:256,:);

    figure('Name','RespVar beha corr topo','Position',[100 20 800 700]);

    for cond_i = 1:3
        for gi = 1:2
            tmpsubs = find(subs.group == gi);
            [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
            subplot(3,3,(gi-1)*3+cond_i);hold all;

            topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
            caxis([-5, 5])
            h = colorbar;
            h.Label.String = sprintf('T-value(p<%.3f)',threshP);
            h.Label.Rotation = -90;
            h.Label.Position = h.Label.Position + [1 0 0];
            title([groupStr{gi},condStr{cond_i}],['NN corr RT TrlCluster=',num2str(Nsamps)]);
        end

        gi =3;% both
        tmpsubs = 1:height(subs);
        [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
        subplot(3,3,(gi-1)*3+cond_i);hold all;

        topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
        caxis([-5, 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i}],['NN corr RT TrlCluster=',num2str(Nsamps)]);

    end

    saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopo_nnSampledAT',num2str(Nsamps),'_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

    %% bar plot of number of significant chans
    myFigBasic;
    myColors = crameri('batlowW',5);
    fig = figure('Position',[300 300 1500 250]);
    subplot(1,4,1);hold all;axis square
    hb = bar(sigChanN.betaChansCorrRT,'FaceAlpha',0,'LineWidth',2);
    for b = 1:length(hb)
        hb(b).FaceColor = myColors(b,:);
        hb(b).EdgeColor = myColors(b,:);
    end
    view(90,90)
    set(gca,'xlim',[0.5 2.5],'XTick',[1 2],'XTickLabel',groupStr)
    ylabel(['Number of channels',newline,'(p<',num2str(threshP),')'],'HorizontalAlignment','center')
    legend(condStr)
    title('Var@NN RT',['TrialCluster=',num2str(Nsamps)])

    % histogram of NN correlation Rho
    myColors = flipud(crameri('bamako',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
    grougAligStr = {'Right','Left'};
    for cond_i = 1:3
        subplot(1,4,cond_i+1);axis square;  hold all;

        for gi  = 1:2
            tmpsubs = find(subs.group == gi);
            [~,p,~,stats] = ttest(subsAll.betaFrontalRTcorr(tmpsubs,cond_i));

            h = histogram(subsAll.betaFrontalRTcorr(tmpsubs,cond_i));
            h.Normalization = 'probability';
            h.BinWidth = 0.15;
            h.FaceColor = myColors(gi,:);
            h.FaceAlpha = 0.9;
            plot([mean(subsAll.betaFrontalRTcorr(tmpsubs,cond_i)) mean(subsAll.betaFrontalRTcorr(tmpsubs,cond_i))],[0 1],'color',myColors(gi,:),'LineWidth',2,'HandleVisibility','off')
            text(2.2-1.1*gi,0.8,sprintf('Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalRTcorr(tmpsubs,cond_i)),stats.tstat,p,computeCohen_d(subsAll.betaFrontalRTcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),...
                'sc','color',myColors(gi,:),'HorizontalAlignment',grougAligStr{gi})
        end
        legend(groupStr,'Location','southeastoutside')
        ylabel('Proportion')
        xlabel('Spearman Rho value')
        title([frontalROI{:}],condStr{cond_i})
        set(gca,'ylim',[0 0.5],'xlim',[-0.6 0.6])
    end
    saveas(gca,fullfile(Dir.figs,['RespBetaVarianceRTchanN_nnSampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
    fig.PaperOrientation = 'landscape';
    print(fig,fullfile(Dir.figs,['RespBetaVarianceRTchanN_nnSampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r300','-bestfit')

    
    %% split young and old based on their NN correlatoin coefficients
    cond_i = 3;
    threshP = 0.05;
    figure('position',[200 200 600 600]);
    for gi = 1:2
        tmpsubs = find(subs.group == gi);
        tmpBeha = subsAll.behaRT(tmpsubs,cond_i);

        tmp_upper = tmpBeha>median(tmpBeha);

        subplot(2,2,(gi-1)*2+1)
        [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs(tmp_upper),:,cond_i)));
        topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
        caxis([-5, 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i},'-highRT'],sprintf('NN corr RT TrlCluster=%d, subjN=%d',Nsamps,sum(tmp_upper)));

        subplot(2,2,(gi-1)*2+2)
        [~,p,~,stats] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs(~tmp_upper),:,cond_i)));
        topoplot(stats.tstat,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
        caxis([-5, 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i},'-lowRT'],sprintf('NN corr RT TrlCluster=%d, subjN=%d',Nsamps,sum(~tmp_upper)));
    end

    saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoSplitRT_nnSampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
end