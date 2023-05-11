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
permN = 10000;
%%
for Nsamps = 40%16:4:40
    clear subsAll
    for sub_i = 1:subN
        subname = subs.name{sub_i};
        matName =fullfile(Dir.results,[subname,'_resp_ft_var_NsampledAT',num2str(Nsamps),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

        if isfile(matName)
            load(matName);

            freqID.beta = dsearchn(tfDat{1}{1}.freq',freq.betaFreq');
            freqID.beta = freqID.beta(1):freqID.beta(2);

            timeID.all = dsearchn(tfDat{1}{1}.time',timeROI.all');
            timeID.all = timeID.all(1):timeID.all(2);

            timeID.Post = dsearchn(tfDat{1}{1}.time',timeROI.Post');
            timeID.Post = timeID.Post(1):timeID.Post(2);

            timeID.Pre = dsearchn(tfDat{1}{1}.time',timeROI.Pre');
            timeID.Pre = timeID.Pre(1):timeID.Pre(2);

            clear chanID
            [log1,chanID.frontal] = ismember(frontalROI,tfDat{1}{1}.label);

            for cond_i = 1:3
                for b = 1:length(tfDat{cond_i})
                    dat{sub_i}.betaAvgFrontal{cond_i}(b,1) = squeeze(mean(mean(mean(tfDat{cond_i}{b}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));

                    dat{sub_i}.betaChans{cond_i}(b,:) = squeeze(mean(mean(tfDat{cond_i}{b}.powspctrm(:,freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"));

                    dat{sub_i}.behaACC{cond_i}(b,1) = mean(PostTrlBeha{cond_i}{b}.button_resp_corr,'omitnan');
                    dat{sub_i}.behaRT{cond_i}(b,1) = mean(PostTrlBeha{cond_i}{b}.button_resp_rt,'omitnan');
                end

                subsAll.betaFrontalVar(sub_i,cond_i) = mean(dat{sub_i}.betaAvgFrontal{cond_i},'omitnan');
                subsAll.behaACC(sub_i,cond_i) = mean(dat{sub_i}.behaACC{cond_i},'omitnan');

                p = polyfit(dat{sub_i}.behaACC{cond_i},dat{sub_i}.betaAvgFrontal{cond_i},1);
                subsAll.betaFrontalACCcorr(sub_i,cond_i) = p(1);

                p = polyfit(dat{sub_i}.behaRT{cond_i},dat{sub_i}.betaAvgFrontal{cond_i},1);
                subsAll.betaFrontalRTcorr(sub_i,cond_i) = p(1);

%                 [subsAll.betaFrontalRTcorr(sub_i,cond_i),subsAll.betaFrontalRT_pval(sub_i,cond_i)]= corr(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.behaRT{cond_i},'type','Spearman');
%                 [subsAll.betaFrontalACCcorr(sub_i,cond_i),subsAll.betaFrontalACC_pval(sub_i,cond_i)] = corr(dat{sub_i}.betaAvgFrontal{cond_i},dat{sub_i}.behaACC{cond_i},'type','Spearman');

                for c = 1:length(tfDat{1}{1}.label)
                    tmp_var = squeeze(dat{sub_i}.betaChans{cond_i}(:,c));
                    tmp_beha = dat{sub_i}.behaACC{cond_i};

                    %                     subsAll.betaChansCorrACC(sub_i,c,cond_i) = corr(tmp_var,tmp_beha,'type','Spearman');
                    p = polyfit(tmp_var,tmp_beha,1);
                    subsAll.betaChansCorrACC(sub_i,c,cond_i) = p(1);

                    tmp_beha =  dat{sub_i}.behaRT{cond_i};
                    %                     subsAll.betaChansCorrRT(sub_i,c,cond_i) = corr(tmp_var,tmp_beha);
                    p = polyfit(tmp_var,tmp_beha,1);
                    subsAll.betaChansCorrRT(sub_i,c,cond_i) = p(1);

                end            
            end           
        end
    end
    save(fullfile(Dir.results,['RespBetaVar_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll','-v7.3')

    %% FrontalBetaACC

% % %     figure('Name','FrontalBetaACC','position',[500 500 900 900]);
% % %     for gi = 1:2
% % %         tmpsubs = find(subs.group == gi);
% % %         [p,t_orig] = mult_comp_perm_t1(subsAll.betaFrontalACCcorr(tmpsubs,:),permN);
% % % 
% % %         for cond_i = 1:3
% % %             subplot(3,3,(gi-1)*3+cond_i);hold all;axis square
% % %             for s = 1:length(tmpsubs)
% % %                 plot(dat{tmpsubs(s)}.betaAvgFrontal{cond_i},dat{tmpsubs(s)}.behaACC{cond_i},'.')
% % %             end
% % %             title(groupStr{gi},[condStr{cond_i},' sampN=' num2str(Nsamps)])
% % %             xlabel(['BetaVariance(' [frontalROI{:}] ')'])
% % %             ylabel('Beha(ACC)')
% % %             set(gca,'ylim',[0 1.5])
% % % 
% % %             text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalACCcorr(tmpsubs,cond_i)),t_orig(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalACCcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
% % %         end
% % %     end
% % % 
% % %     gi =3;% both groups
% % %     tmpsubs = 1:height(subs);
% % %     [p,t_orig] = mult_comp_perm_t1(subsAll.betaFrontalACCcorr(tmpsubs,:),permN);
% % %     for cond_i = 1:3
% % %         subplot(3,3,(gi-1)*3+cond_i);hold all;axis square
% % %         for s = 1:length(tmpsubs)
% % %             plot(dat{tmpsubs(s)}.betaAvgFrontal{cond_i},dat{tmpsubs(s)}.behaACC{cond_i},'.')
% % %             lsline
% % %         end
% % %         title(groupStr{gi},[condStr{cond_i},' sampN=' num2str(Nsamps)])
% % %         xlabel(['BetaVariance(' [frontalROI{:}] ')'])
% % %         ylabel('Beha(ACC)')
% % %         set(gca,'ylim',[0 1.5])
% % % 
% % %         text(0.95,0.85,sprintf('mean Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalACCcorr(tmpsubs,cond_i)),t_orig(cond_i),p(cond_i),computeCohen_d(subsAll.betaFrontalACCcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),'sc','HorizontalAlignment','Right')
% % %     end
% % %     saveas(gca,fullfile(Dir.figs,['RespBetaVarianceFrontalACC_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'_perm.png']))

      %% correlation topography beta

    myFigBasic;

    load chanLocs_dms.mat
    [is, loc] = ismember(tfDat{1}{1}.label,{chanloc.labels});
    chanloc  = chanloc(loc(is));
    threshP = 0.05;
    mkSize = 8;
    myColors = jet;
%     myColors = hot;
    %     myColors = myColors(256/2:256,:);
    % myColors = crameri('bam',256);
    figure('Name','RespVar beha corr topo','Position',[100 20 800 1500]);

    for cond_i = 1:3
        for gi = 1:2
            tmpsubs = find(subs.group == gi);
%             [p,t_orig] = mult_comp_perm_t1(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)),permN);
           [~,p,~,stat] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
t_orig = stat.tstat;

subplot(6,3,(gi-1)*3+cond_i);hold all;

            topoplot(t_orig,chanloc,'emarker2',{find(p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
            caxis([-5, 5])
            h = colorbar;
            h.Label.String = sprintf('T-value(p<%.3f)',threshP);
            h.Label.Rotation = -90;
            h.Label.Position = h.Label.Position + [1 0 0];
            title([groupStr{gi},condStr{cond_i}],['N+1 corr RT TrlCluster=',num2str(Nsamps)]);

%             [p,t_orig] = mult_comp_perm_t1(squeeze(subsAll.betaChansCorrACC(tmpsubs,:,cond_i)),permN);
            
[~,p,~,stat] = ttest(squeeze(subsAll.betaChansCorrACC(tmpsubs,:,cond_i)));
t_orig = stat.tstat;

            sigChanN.betaChansCorrACC(gi,cond_i) = sum((p<threshP & t_orig>0));% for bar plot
            sigChanN.betaChansCorrACC(gi,cond_i) = sum(p<threshP);% for bar plot

            subplot(6,3,(gi-1)*3+cond_i+9);hold all;

            topoplot(t_orig,chanloc,'emarker2',{find(p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
            caxis([-3, 3])
            h = colorbar;
            h.Label.String = sprintf('T-value(p<%.3f)',threshP);
            h.Label.Rotation = -90;
            h.Label.Position = h.Label.Position + [1 0 0];
            title([groupStr{gi},condStr{cond_i}],['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
        end

        gi =3;% both
        tmpsubs = 1:height(subs);
%         [p,t_orig] = mult_comp_perm_t1(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)),permN);

[~,p,~,stat] = ttest(squeeze(subsAll.betaChansCorrRT(tmpsubs,:,cond_i)));
t_orig = stat.tstat;

        subplot(6,3,(gi-1)*3+cond_i);hold all;

        topoplot(t_orig,chanloc,'emarker2',{find(p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
        caxis([-5, 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i}],['N+1 corr RT TrlCluster=',num2str(Nsamps)]);

%         [p,t_orig] = mult_comp_perm_t1(squeeze(subsAll.betaChansCorrACC(tmpsubs,:,cond_i)),permN);

[~,p,~,stat] = ttest(squeeze(subsAll.betaChansCorrACC(tmpsubs,:,cond_i)));
t_orig = stat.tstat;

        subplot(6,3,(gi-1)*3+cond_i+9);hold all;

        topoplot(t_orig,chanloc,'emarker2',{find(p<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
        caxis([-5, 5])
        h = colorbar;
        h.Label.String = sprintf('T-value(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title([groupStr{gi},condStr{cond_i}],['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
    end

    saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopo_SampledAT',num2str(Nsamps),'_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '_perm.png']))

    %% 2*3 mixed ANOVA on correlation coefficients
%     threshP = .005;
%     clear xFs xPs
%     for c = 1:size(subsAll.betaChansCorrACC,2)
%         tmpdat = squeeze(subsAll.betaChansCorrACC(:,c,:)); % sub*chan*load
%         tmpdat = reshape(tmpdat,size(tmpdat,1)*size(tmpdat,2),1);% [sub*load]*1
% 
%         clear X
%         %     X(:,1) = tmpdat;
%         X(:,1) = atanh(tmpdat); % fisherZ transformation of correlation coefficients
%         X(:,2) = flipud(repmat(subs.group,3,1));
%         X(:,[4 3]) = fullfact([subN 3]);
% 
%         [~, ~, ~, xFs{c}, xPs{c}]=mixed_between_within_anova(X);
%     end
% 
%     clear xFsMat
%     for c = 1:size(subsAll.betaChansCorrACC,2)
%         xFsMat.inter(c) = xFs{c}{4};
%         xFsMat.load(c) = xFs{c}{3};
%         xFsMat.group(c) = xFs{c}{1};
% 
%         xFsMat.inter_p(c) = xPs{c}{4};
%         xFsMat.load_p(c) = xPs{c}{3};
%         xFsMat.group_p(c) = xPs{c}{1};
% 
%     end
%     figure('Position',[200 200 800 300]);
%     for i = 1:3
%         subplot(1,3,i);
%         if i ==1
%             tmpdat = xFsMat.inter;
%             tmp_p = xFsMat.inter_p;
%             title('Interaction effect',['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
%         elseif  i==2
%             tmpdat = xFsMat.load;
%             tmp_p = xFsMat.load_p;
%             title('Load effect effect',['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
% 
%         elseif i==3
%             tmpdat = xFsMat.group;
%             tmp_p = xFsMat.group_p;
%             title('Group effect effect',['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
% 
%         end
%         topoplot(tmpdat,chanloc,'emarker2',{find(tmp_p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
%         caxis([0, 5])
%         h = colorbar;
%         h.Label.String = sprintf('T-value(p<%.3f)',threshP);
%         h.Label.Rotation = -90;
%         h.Label.Position = h.Label.Position + [1 0 0];
%     end
%     saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoANOVA_SampledAT',num2str(Nsamps),'_threshP',num2str(threshP),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '_perm.png']))
%     
    %% bar plot of number of significant chans
    myFigBasic;
    myColors = crameri('batlowW',5);
    fig = figure('Position',[300 300 1500 250]);
    subplot(1,4,1);hold all;axis square
    hb = bar(sigChanN.betaChansCorrACC,'FaceAlpha',0,'LineWidth',2);
    for b = 1:length(hb)
        hb(b).FaceColor = myColors(b,:);
        hb(b).EdgeColor = myColors(b,:);
    end
    view(90,90)
    set(gca,'xlim',[0.5 2.5],'XTick',[1 2],'XTickLabel',groupStr)
    ylabel(['Number of channels',newline,'(p<',num2str(threshP),')'],'HorizontalAlignment','center')
    legend(condStr)
    title('Var@N+1 ACC',['TrialCluster=',num2str(Nsamps)])

    % histogram of N+1 correlation Rho
    myColors = flipud(crameri('bamako',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
    grougAligStr = {'Right','Left'};
    for cond_i = 1:3
        subplot(1,4,cond_i+1);axis square;  hold all;

        for gi  = 1:2
            tmpsubs = find(subs.group == gi);
            [p,t_orig] = mult_comp_perm_t1(subsAll.betaFrontalACCcorr(tmpsubs,cond_i),permN);

            h = histogram(subsAll.betaFrontalACCcorr(tmpsubs,cond_i));
            h.Normalization = 'probability';
            h.BinWidth = 0.15;
            h.FaceColor = myColors(gi,:);
            h.FaceAlpha = 0.9;
            plot([mean(subsAll.betaFrontalACCcorr(tmpsubs,cond_i)) mean(subsAll.betaFrontalACCcorr(tmpsubs,cond_i))],[0 1],'color',myColors(gi,:),'LineWidth',2,'HandleVisibility','off')
            text(2.2-1.1*gi,0.8,sprintf('Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.betaFrontalACCcorr(tmpsubs,cond_i)),t_orig,p,computeCohen_d(subsAll.betaFrontalACCcorr(tmpsubs,cond_i),zeros(length(tmpsubs),1),'paired')),...
                'sc','color',myColors(gi,:),'HorizontalAlignment',grougAligStr{gi})
        end
        legend(groupStr,'Location','southeastoutside')
        ylabel('Proportion')
        xlabel('Spearman Rho value')
        title([frontalROI{:}],condStr{cond_i})
        set(gca,'ylim',[0 0.5],'xlim',[-0.6 0.6])
    end
    saveas(gca,fullfile(Dir.figs,['RespBetaVarianceACCchanN_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'_perm.png']))
    fig.PaperOrientation = 'landscape';
    print(fig,fullfile(Dir.figs,['RespBetaVarianceACCchanN_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r300','-bestfit')
 
    %% split young and old based on their N+1 correlatoin coefficients
%     cond_i = 3;
%     threshP = 0.05;
%     figure('position',[200 200 600 600]);
%     for gi = 1:2
%         tmpsubs = find(subs.group == gi);
%         tmpBeha = subsAll.behaACC(tmpsubs,cond_i);
% 
%         tmp_upper = tmpBeha>median(tmpBeha);
% 
%         subplot(2,2,(gi-1)*2+1)
%         [p,t_orig] = mult_comp_perm_t1(squeeze(subsAll.betaChansCorrACC(tmpsubs(tmp_upper),:,cond_i)),permN);
%         topoplot(t_orig,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
%         caxis([-5, 5])
%         h = colorbar;
%         h.Label.String = sprintf('T-value(p<%.3f)',threshP);
%         h.Label.Rotation = -90;
%         h.Label.Position = h.Label.Position + [1 0 0];
%         title([groupStr{gi},condStr{cond_i},'-highACC'],sprintf('N+1 corr ACC TrlCluster=%d, subjN=%d',Nsamps,sum(tmp_upper)));
% 
%         subplot(2,2,(gi-1)*2+2)
%         [p,t_orig] = mult_comp_perm_t1(squeeze(subsAll.betaChansCorrACC(tmpsubs(~tmp_upper),:,cond_i)),permN);
%         topoplot(t_orig,chanloc,'emarker2',{find(p<threshP),'p','k',mkSize},'plotrad',0.53,'electrodes','off' );
%         caxis([-5, 5])
%         h = colorbar;
%         h.Label.String = sprintf('T-value(p<%.3f)',threshP);
%         h.Label.Rotation = -90;
%         h.Label.Position = h.Label.Position + [1 0 0];
%         title([groupStr{gi},condStr{cond_i},'-lowACC'],sprintf('N+1 corr ACC TrlCluster=%d, subjN=%d',Nsamps,sum(~tmp_upper)));
%     end
% 
%     saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoSplitACC_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '_perm.png']))

    %% spatial cluster-permutation on N+1 correlation
    threshP = 0.05;
    for gi = 1:2
        ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));
        cfg= [];
        cfg.latency = 0;
        cfg.keepindividual = 'yes';
        tmpsubs = subs.group==gi;
        cfg.channel = tfDat{1}{1}.label;
        ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{tmpsubs,1});

        cfg = [];
        cfg.layout = 'easycapM11.mat';
        layout = ft_prepare_layout(cfg);

        cfg = [];
        cfg.layout = layout;
        cfg.method = 'triangulation';
        % cfg.compress = 'yes';
        % cfg.feedback = 'yes';
        neighbours = ft_prepare_neighbours(cfg);

        for cond_i = 1:3
            ft_exp = ft_temp;
            ft_exp.individual = squeeze(subsAll.betaChansCorrACC(tmpsubs,:,cond_i));

            ft_zero = ft_temp;
            ft_zero.individual = ft_zero.individual*0;

            cfg                  = [];
            %         cfg.channel = 'All';
            cfg.parameter = 'individual';
            %         cfg_neighb.method = 'distance';
            cfg.neighbours       = neighbours;
            cfg.method           = 'montecarlo';
            cfg.statistic        = 'indepsamplesT';
            cfg.correctm         = 'cluster';

            cfg.clusteralpha     = threshP;
            cfg.clusterstatistic = 'maxsum';
            cfg.minnbchan        = 2;%% minimal number of neighbouring channels
            cfg.tail             = 0;
            cfg.clustertail      = 0;
            cfg.alpha            = threshP;
            cfg.numrandomization = 10000;

            cfg.design           = [ones(1,sum(subs.group==gi)), ones(1,sum(subs.group==gi))*2]; % design matrix
            cfg.ivar             = 1; % number or list with indices indicating the independent variable(s)
            chanStat{gi,cond_i}  = ft_timelockstatistics(cfg, ft_exp, ft_zero);
        end
    end

    %%  make a plot
    myColors = hot;
    for gi = 1:2
        for cond_i = 1:3
            figure('position',[100 100 250 250]);
            cfg = [];
            cfg.layout = 'easycapM11.mat';
            cfg.parameter='stat';
            cfg.contournum = 0;
            cfg.zlim = [-6 6];
            cfg.comment = 'no';
            if sum(chanStat{gi,cond_i}.mask)>0 &&  ~isempty(chanStat{gi,cond_i}.posclusters)

                cfg.maskparameter = 'mask';
                cfg.highlightsymbolseries = ['.','.','.','.','.'];
                cfg.highlightsymbol  = '.';
                cfg.highlightsize=20;
                %                 cfg.markersymbol = '.';
                cfg.subplotsize = [1 1];
                ft_clusterplot(cfg, chanStat{gi,cond_i} );
            else
                cfg.markersize=0.000001;
                ft_topoplotER(cfg, chanStat{gi,cond_i})
            end

            colormap(myColors);
            h = colorbar;
            h.Label.String = sprintf('T-value(cluster permutation with p<%.3f)',threshP);
            h.Label.Rotation = -90;
            h.Label.Position = h.Label.Position + [1 0 0];
            title(groupStr{gi},[condStr{cond_i},' sampN=' num2str(Nsamps)])

            saveas(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoN+1permutation_SampledAT',num2str(Nsamps),groupStr{gi},condStr{cond_i},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
            print(gcf,fullfile(Dir.figs,['RespBetaVarianceTopoN+1permutation_SampledAT',num2str(Nsamps),groupStr{gi},condStr{cond_i},num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']),'-dpdf','-r300','-bestfit')
        end
    end

    myColors = crameri('batlowW',5);
    fig = figure('Position',[300 300 250 250]);
    hold all;axis square

    tmp = @(x)sum(x.mask);
    chanN = cellfun(tmp,chanStat);

    hb = bar(chanN,'FaceAlpha',1,'LineWidth',0.5);
    for b = 1:length(hb)
        hb(b).FaceColor = myColors(b,:);
%         hb(b).EdgeColor = myColors(b,:);
    end
    view(90,90)
    set(gca,'xlim',[0.5 2.5],'XTick',[1 2],'XTickLabel',groupStr)
    ylabel(['Number of channels',newline,'(cluster permutation p<',num2str(threshP),')'],'HorizontalAlignment','center')
    legend(condStr)
    title('Var@N+1 ACC',['TrialCluster=',num2str(Nsamps)])

    saveas(gca,fullfile(Dir.figs,['RespBetaVarianceACCchanNpermed_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'_perm.png']))
    saveas(gca,fullfile(Dir.figs,['RespBetaVarianceACCchanNpermed_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'_perm.pdf']))
end