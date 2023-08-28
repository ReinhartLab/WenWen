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
IsBL2preprobe = 0; % use pre-trial baseline
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preprobe'};
groupStr = {'Young','Old','Young-Old'};
condStr = {'ss1','ss2','ss4'};
condDiffStr = {'S2-S1','S4-S1','S4-S2'};

frontalROI = {'Fz','F1','F2','FCz','FC1','FC2'};
% frontalROI = {'Fz','F1','F2',};
occipROI = {'Pz','POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz
timeROI = [-1 3];% in s
SStime = {[-1.6 4.5],[-2.8 4.5],[-5.2 4.5]};% epoch length for different SS

%%
subsAll = cell(subN,3);
clear tmpConn
for sub_i = 1:subN
    subname = subs.name{sub_i};

    matName =fullfile(Dir.results,[subname,'_probe_ft_plv',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preprobe+1,4},'.mat']);

    if isfile(matName)
        load(matName);
        for ss = 1:3
            %             tfDat.plv{ss}.plvspctrm = bsxfun(@rdivide,tfDat.plv{ss}.plvspctrm, tfDat.plvB{ss}.plvspctrm);
            tmpConn(sub_i,ss,:,:,:,:) = tfDat.plv{ss}.plvspctrm;
        end

        subsAll(sub_i,:) = tfDat.plv;
    end
end

empty_idx = cellfun(@isempty,subsAll(:,1));
subsAll(empty_idx,:) = [];
subs(empty_idx,:) = [];
%%
for gi = 1:2
    snList = subs.group == gi;
    for ss = 1:3
        gndTF{gi,ss} = subsAll{sub_i,ss};% get structure
        gndTF{gi,ss}.plvspctrm =  squeeze(mean(tmpConn(snList,ss,:,:,:,:),'omitnan'));
    end
end

for gi = 1:2
    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'plvspctrm';
    gndTF_SSdiff{gi,1} = ft_math(cfg, gndTF{gi,2},gndTF{gi,1}); %set size difference
    gndTF_SSdiff{gi,2} = ft_math(cfg, gndTF{gi,3},gndTF{gi,1});
    gndTF_SSdiff{gi,3} = ft_math(cfg, gndTF{gi,3},gndTF{gi,2});
end

for ss = 1:3
    cfg = [];
    cfg.operation='subtract';
    cfg.parameter = 'plvspctrm';
    gndTF_Groupdiff{ss} = ft_math(cfg, gndTF{1,ss},gndTF{2,ss}); %group difference
end

%% time-frequency map
myFigBasic;

w = [1 1 1];
b = [0,0,1];
k = [0,0,0];
r = [1,0,0];
bkr = createcolormap(w,b,k,r,w); % 256x3 array

cfg = [];
cfg.xlim  = [-0.4 0.7];
% cfg.ylim = [1 40];
cfg.zlim = [0 0.25];
cfg.parameter = 'plvspctrm';
cfg.colormap = turbo;
cfg.refchannel = {'T7'};
cfg.channel = {'Fz','F1','F3','FC1'};
cfg.figure = 'gca';
cfg.colorbar  = 'no';

figure('Position',[100 100 1400 600]);
for gi = 1:2
    for cond_i = 1:3
        subplot(2,3,(gi-1)*3+cond_i);hold all;axis square
        ft_singleplotTFR(cfg, gndTF{gi,cond_i});
% colorbar;
% hc = colorbar;
%         hc.Label.String = 'PLV';
%         hc.Label.Rotation = -90;
%         hc.Label.Position = [4.2 0 0];

        title([condStr{cond_i},groupStr{gi}],[cfg.channel{:}],'FontSize',14);colorbar
ylabel('\bfFrequency(Hz)');
        xlabel('Time(0s=Probe)')
            rectangle('Position',[0 0 0 40],'LineWidth',1,'LineStyle',':')
end
end
   
       
% saveas(gcf,fullfile(Dir.figs,['PLVprobe_',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preprobe+1,4} '.bmp']))


%% time-resolved topopplot based on ref channel
load chanLocs_dms.mat
chanlocs = chanloc(ismember({chanloc.labels},gndTF{1}.elec.label));
cfg = [];
cfg.ylim = [3 7];
cfg.refchannel = 'T7';
cfg.zlim = [-0.02 0.02];

cfg.timeBin = {[0 0.1],[0.1 0.2],[0.2 0.3],[0.3 0.4],[0.4 0.5],[0 0.5]};
figure('Position',[10 10 2000 1200])
for ti = 1:length(cfg.timeBin)
    cfg.xlim = cfg.timeBin{ti};

    for gi = 1:2
        [~,refchanID] = ismember(cfg.refchannel,gndTF{gi,1}.label);
        tmpfreq = dsearchn(gndTF{gi,1}.freq',cfg.ylim');
        tmptoi = dsearchn(gndTF{gi,1}.time',cfg.xlim');

        tmp = diff(tmpConn{gi}(:,[1 2],:,:,:,:),1,2);
        subdatTopo{gi}(1,:,:) = squeeze(mean(mean(tmp(:,1,refchanID,:,tmpfreq,tmptoi),5),6));

        tmp = diff(tmpConn{gi}(:,[1 3],:,:,:,:),1,2);
        subdatTopo{gi}(2,:,:) = squeeze(mean(mean(tmp(:,1,refchanID,:,tmpfreq,tmptoi),5),6));

        tmp = diff(tmpConn{gi}(:,[2 3],:,:,:,:),1,2);
        subdatTopo{gi}(3,:,:) = squeeze(mean(mean(tmp(:,1,refchanID,:,tmpfreq,tmptoi),5),6));

        for cond_i = 1:3
            subplot(length(cfg.timeBin),6,(ti-1)*6+cond_i+(gi-1)*3);axis square

            hmask{gi,cond_i} = ttest(squeeze(subdatTopo{gi}(cond_i,:,:)),0,'alpha',0.01);

            topoplot(squeeze(mean(subdatTopo{gi}(cond_i,:,:),2)),chanlocs,'plotrad',.53,'electrodes','off','style','map','emarker2',{find(hmask{gi,cond_i}),'p','k'},'maplimits',[-.02 .02]);
            colormap(bkr);
            hc = colorbar;
            hc.Label.String = 'PLV';
            hc.Label.Rotation = -90;
            hc.Label.Potation = [4 0 0  ];
            title(sprintf('%s\n%s',groupStr{gi},condDiffStr{cond_i}),[num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)) 'Hz: ' num2str(cfg.xlim(1)),'~',num2str(cfg.xlim(2)) 's']);
        end
    end
end
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_topo_',num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)),'Hz_ref',cfg.refchannel,txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_topo_',num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)),'Hz_ref',cfg.refchannel,txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.emf']))
%%
figure('Position',[10 10 900 450])
ti = length(cfg.timeBin);
cfg.xlim = cfg.timeBin{ti};

for gi = 1:2
    [~,refchanID] = ismember(cfg.refchannel,gndTF{gi,1}.label);
    tmpfreq = dsearchn(gndTF{gi,1}.freq',cfg.ylim');
    tmptoi = dsearchn(gndTF{gi,1}.time',cfg.xlim');

    tmp = diff(tmpConn{gi}(:,[1 2],:,:,:,:),1,2);
    subdatTopo{gi}(1,:,:) = squeeze(mean(mean(tmp(:,1,refchanID,:,tmpfreq,tmptoi),5),6));

    tmp = diff(tmpConn{gi}(:,[1 3],:,:,:,:),1,2);
    subdatTopo{gi}(2,:,:) = squeeze(mean(mean(tmp(:,1,refchanID,:,tmpfreq,tmptoi),5),6));

    tmp = diff(tmpConn{gi}(:,[2 3],:,:,:,:),1,2);
    subdatTopo{gi}(3,:,:) = squeeze(mean(mean(tmp(:,1,refchanID,:,tmpfreq,tmptoi),5),6));

    for cond_i = 1:3
        subplot(2,3,cond_i+(gi-1)*3);axis square

        hmask{gi,cond_i} = ttest(squeeze(subdatTopo{gi}(cond_i,:,:)),0,'alpha',0.05);

        topoplot(squeeze(mean(subdatTopo{gi}(cond_i,:,:),2)),chanlocs,'plotrad',.53,'electrodes','off','style','map','emarker2',{find(hmask{gi,cond_i}),'p','w',8},'maplimits',[-.022 .022]);
        colormap(bkr);
        hc=colorbar;
        hc.FontSize = 14;
        hc.Label.String = 'PLV';
        hc.Label.Position = [5,0,0];
        hc.Label.Rotation = -90;

        %         title(sprintf('%s\n%s',groupStr{gi},condDiffStr{cond_i}),[num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)) 'Hz: ' num2str(cfg.xlim(1)),'~',num2str(cfg.xlim(2)) 's'],'FontSize',14);
        title(sprintf('%s\n%s',groupStr{gi},condDiffStr{cond_i}),'FontSize',14);
    end
end

% saveas(gcf,fullfile('Figs',['ProbePLV_ft_topo_avg',num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)),'Hz_ref',cfg.refchannel,txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_topo_avg',num2str(cfg.ylim(1)),'-',num2str(cfg.ylim(2)),'Hz_ref',cfg.refchannel,txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.pdf']))


%% correlation
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bam',4);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

addpath(genpath('D:\Toolbox\gramm-master'))
load Behav_Y.mat
wanted = ismember(subs.name(mySubs{1}),demoInfo.name);
% dat.beha{1} = squeeze(B.swapNT(wanted,:,1));% k of fixation ignore and attend condition
% dat.beha{1} = squeeze(B.swapDots(wanted,:,1));
% dat.beha{1} = memResults.rotBegin(wanted,:);
% dat.beha{1} = memResults.rotACC(wanted,:);
dat.beha{1} = memResults.rotDif(wanted,:);

load Behav_E.mat
wanted = ismember(subs.name(mySubs{2}),demoInfo.name);
% dat.beha{2} = squeeze(B.swapNT(wanted,:,1));% k of fixation ignore and attend condition
% dat.beha{2} = squeeze(B.swapDots(wanted,:,1));
% dat.beha{2} = memResults.rotBegin(wanted,:);
% dat.beha{2} = memResults.rotACC(wanted,:);
dat.beha{2} = memResults.rotDif(wanted,:);

% load DecodingProbeTacc.mat

clear g
for gi = 1:2
    [r,pval] = corr(ConnAvg{gi},dat.beha{gi},'type','Spearman','rows','complete');

    for cond = 1:4
        Y = dat.beha{gi}(:,cond);
        X = ConnAvg{gi}(:,cond);

        g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,cond).geom_point(); %Draw a scatter plot
        g(gi,cond).stat_glm('disp_fit',0); %Draw lines and confidence intervals based on scatter plots
        g(gi,cond).set_names('x','PLV','y','Recall difference'); %Set the title of the axis
        g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
        g(gi,cond).set_title(sprintf('%s: Rho = %.3f, p = %.3f',groupStr{gi},r(cond,cond),pval(cond,cond)));
        g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    end
end
figure('Position',[100 100 1000 500]);
g.axe_property('XLim',[0 0.25]);
g.draw();
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_gramm_',num2str(avgFreq(1)),'-',num2str(avgFreq(2)),'Hz_ref',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2)),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_gramm_',num2str(avgFreq(1)),'-',num2str(avgFreq(2)),'Hz_ref',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2)),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.emf']))

%% frequency-resolved 
myColors = crameri('devon',4);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

clear freqPeak
refchan = 'T7';
refchanID = find(strcmp(refchan,gndTF{1}.label));
tochan = {'F1','F3','Fz','FC1'};
[~,tochanID] = ismember(tochan,gndTF{1}.label);

avgToi = [0 0.5];
tmptoi = dsearchn(gndTF{1}.time',avgToi');

for gi = 1:2
    snList = subs.group == gi;
    freqPeak{gi} = squeeze(mean(mean(tmpConn(snList,:,refchanID,tochanID,:,tmptoi),4,'omitnan'),6,'omitnan'));
end
% interaction effect of mixed anova
freqPeak_pval = nan(size(freqPeak{1},3),1);
for fi = 1:size(freqPeak{1},3)
    
    if isnan(freqPeak{1}(1,1,fi))
        continue
    end
    X = [];
    for gi = 1:2
        clear tmp
        tmp(:,1) = reshape(freqPeak{gi}(:,:,fi),size(freqPeak{gi},1)*size(freqPeak{gi},2),[]);
        tmp(:,2) = ones(length(tmp),1)*gi;
        tmp(:,3) = sort(repmat(1:size(freqPeak{gi},2),1,size(freqPeak{gi},1))');
        tmp(:,4) = repmat(1:size(freqPeak{gi},1),1,size(freqPeak{gi},2))'+(gi-1)*size(freqPeak{1},1);
        X = [X;tmp];
    end

 [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);

freqPeak_pval(fi) = Ps{4};
end

freqX = gndTF{1}.freq;

clear PlotSubj
% PlotSubj.name = 'SE19';
% [~,PlotSubj.ID] = ismember(PlotSubj.name,subs.info.name);
% PlotSubj.group = subs.info.groupStr(PlotSubj.ID);

figure('Position',[200 200 750 350]);
linS = {'-.','-'};

for gi = 1:2
    subplot(1,3,gi);axis square;hold all;box on

    mn = squeeze(mean(freqPeak{gi},'omitnan'))';
    se = squeeze(std(freqPeak{gi},'omitnan'))'./sqrt(sum(subs.group == gi));

    for i = 1:size(mn,2)
%         hl(i).Color = myColors(i,:);
        shadedErrorBar(freqX,mn(:,i),se(:,i),{'color',myColors(i,:),'LineStyle',linS{gi},'linewidth',1.2})
    end

    if exist('PlotSubj','var')
        ismember(PlotSubj.ID,mySubs{gi})
        tmpDat = squeeze(freqPeak{gi}(mySubs{gi}==PlotSubj.ID,:,:));
        hl = plot(tmpDat','--');

        for i = 1:length(hl)
            hl(i).Color = myColors(i,:);
        end
    end

    set(gca,'YLim',[0.1 0.3],'xTick',[0:5:40])

    % rectangle('Position',[0 0 10 0.5],'FaceColor',[0.5 0.5 0.5 0.2],'EdgeColor',[1 1 1 0])
    ytickformat('%.2f')

    legend(condStr)
    title(groupStr{gi})
    xlabel('Frequency(Hz)')
    ylabel('PLV')
end

subplot(1,3,3);axis square;hold all;box on
for gi = 1:2
    mn = squeeze(mean(freqPeak{gi},'omitnan'))';
     hl = plot(mn,'LineStyle',linS{gi},'LineWidth',1.2);
     for i = 1:length(hl)
         hl(i).Color = myColors(i,:);
%          shadedErrorBar(freqX,mn(:,i),se(:,i),{'color',myColors(i,:),'LineStyle',linS{gi}})
    end
end
%plot significance
text(find(freqPeak_pval<=.05),0.09*ones(sum(freqPeak_pval<=.05),1),'*','HorizontalAlignment','center','FontSize',8)

ytickformat('%.2f')
tmp = dsearchn(gndTF{1,1}.freq',[0:5:40]');
set(gca,'YLim',[0.1 0.3],'YTick',[0.08:0.04:0.2],'xTick',tmp,'xTicklabel',[0:5:40])
tmp = dsearchn(gndTF{1,1}.freq',[3 7]');
rectangle('Position',[tmp(1) 0 diff(tmp) 0.5],'FaceColor',[0.5 0.5 0.5 0.2],'EdgeColor',[1 1 1 0])
title('Both')
xlabel('Frequency(Hz)')
ylabel('PLV')

% % axes('Position',[0.68 0.49 0.32 0.24],'Box','on')
% % hold all;axis square
% % dat2plot = cellfun(@mean,ConnAvg,'UniformOutput',false);
% % dat2plot = vertcat(dat2plot{:});%group*cond
% % 
% % se2plot = cellfun(@std,ConnAvg,'UniformOutput',false);
% % se2plot = vertcat(se2plot{:})./sqrt(cellfun(@length,ConnAvg))';%group*cond
% % 
% % hb = bar(dat2plot);
% % for i = 1:length(hb)
% %     hb(i).FaceColor = myColors(i,:);
% % end
% % xcord = vertcat(hb.XEndPoints)';
% % errorbar(xcord,dat2plot,se2plot,'k','LineStyle','none')
% % set(gca,'xtick',[1 2],'xticklabel',groupStr,'YLim',[0.1 0.2],'YTick',[0.1 0.2]);
% title([num2str(avgFreq(1)),'-',num2str(avgFreq(2)),'Hz_ref',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2))],'Interpreter','none')

% [~,p] = ttest(diff(ConnAvg{1}(:,[2 3]),1,2))
% plot(xcord(1,[1 3]),[1 1]*0.16,'k','HandleVisibility','off')
% text(mean(xcord(1,[1 3])),0.16,'**','FontSize',14,'HorizontalAlignment','center');
% plot(xcord(1,[2 3]),[1 1]*0.17,'k','HandleVisibility','off')
% text(mean(xcord(1,[2 3])),0.17,'***','FontSize',14,'HorizontalAlignment','center');
% 
% plot(xcord(1,[1 4]),[1 1]*0.18,'k','HandleVisibility','off')
% text(mean(xcord(1,[1 4])),0.18,'***','FontSize',14,'HorizontalAlignment','center');
% plot(xcord(1,[2 4]),[1 1]*0.19,'k','HandleVisibility','off')
% text(mean(xcord(1,[2 4])),0.19,'***','FontSize',14,'HorizontalAlignment','center');

% saveas(gcf,fullfile('Figs',['ProbePLV_ft_freqPeak_ref',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2)),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_freqPeak_ref',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2)),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.emf']))
% %% for SFN poster
% myFigBasic
% figure('Position',[100 100 790 300])
% subplot(121);hold all;box on;axis square
% for gi = 1:2
%     mn = squeeze(mean(freqPeak{gi}))';
%     plot(mn(:,1),'LineStyle',linS{gi},'color','k')
%     hl = plot(mn,'LineStyle',linS{gi},'HandleVisibility','off');
%     for i = 1:length(hl)
%         hl(i).Color = myColors(i,:);
%         %          shadedErrorBar(freqX,mn(:,i),se(:,i),{'color',myColors(i,:),'LineStyle',linS{gi}})
%     end
% end
% text(find(freqPeak_pval<=.05),0.09*ones(sum(freqPeak_pval<=.05),1),'*','HorizontalAlignment','center','FontSize',12)
% 
% ytickformat('%.2f')
% hl = legend(groupStr);
% tmp = dsearchn(gndTF{1,1}.freq',[0:5:40]');
% set(gca,'YLim',[0.08 0.181],'ytick',[0.1:0.02:0.2],'xTick',tmp,'xTicklabel',[0:5:40])
% tmp = dsearchn(gndTF{1,1}.freq',[3 7]');
% rectangle('Position',[tmp(1) 0 diff(tmp) 0.5],'FaceColor',[0.5 0.5 0.5 0.4],'EdgeColor',[1 1 1 0])
% xlabel('Frequency(Hz)')
% ylabel('PLV')
% 
% subplot(122);
% hold all;box on
% dat2plot = cellfun(@mean,ConnAvg,'UniformOutput',false);
% dat2plot = vertcat(dat2plot{:});%group*cond
% 
% se2plot = cellfun(@std,ConnAvg,'UniformOutput',false);
% se2plot = vertcat(se2plot{:})./sqrt(cellfun(@length,ConnAvg))';%group*cond
% 
% hb = bar(dat2plot);
% 
% for i = 1:length(hb)
%     hb(i).FaceColor = myColors(i,:);
% end
% xcord = vertcat(hb.XEndPoints)';
% errorbar(xcord,dat2plot,se2plot,'k','LineStyle','none')
% set(gca,'xtick',[1 2],'xticklabel',groupStr,'XLim',[0.5 2.5],'YLim',[0.1 0.2],'YTick',[0.1:0.02: 0.2]);
% % title([num2str(avgFreq(1)),'-',num2str(avgFreq(2)),'Hz_ref',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2))],'Interpreter','none')
% ytickformat('%.2f')
% tmp = linspace(0.165,0.185,4);
% plot(xcord(1,[1 3]),[1 1]*tmp(1),'k','HandleVisibility','off')
% text(mean(xcord(1,[1 3])),tmp(1),'**','FontSize',20,'HorizontalAlignment','center');
% plot(xcord(1,[2 3]),[1 1]*tmp(2),'k','HandleVisibility','off')
% text(mean(xcord(1,[2 3])),tmp(2),'***','FontSize',20,'HorizontalAlignment','center');
% 
% plot(xcord(1,[1 4]),[1 1]*tmp(3),'k','HandleVisibility','off')
% text(mean(xcord(1,[1 4])),tmp(3),'***','FontSize',20,'HorizontalAlignment','center');
% plot(xcord(1,[2 4]),[1 1]*tmp(4),'k','HandleVisibility','off')
% text(mean(xcord(1,[2 4])),tmp(4),'***','FontSize',20,'HorizontalAlignment','center');
% legend(condStr)
% saveas(gcf,fullfile('Figs',['ProbePLV_ft_freqPeak_ref_sfn',refchan,'_to_',[tochan{:}],'_',num2str(avgToi(1)),'~',num2str(avgToi(2)),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.emf']))


%%

freqID.beta = dsearchn(gndTF{1}.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

freqID.alpha = dsearchn(gndTF{1}.freq',freq.alphaFreq');
freqID.alpha = freqID.alpha(1):freqID.alpha(2);

timeID = dsearchn(gndTF{1}.time',timeROI');
timeID = timeID(1):timeID(2);

clear chanID
[log1,chanID.frontal] = ismember(frontalROI,gndTF{1}.label);
[log2,chanID.occip] = ismember(occipROI,gndTF{1}.label);

clear dat

for gi = 1:2
    for cond_i = 1:3
        dat.betaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID),2),3),4));
        dat.betaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID),2),3),4));

        dat.alphaAvgOccip{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID),2),3),4));
        dat.alphaAvgFrontal{gi}(:,cond_i) = squeeze(mean(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID),2),3),4));

        dat.betaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.beta,timeID),2),3));
        dat.betaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.beta,timeID),2),3));

        dat.alphaCurveOccip{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.occip(log2),freqID.alpha,timeID),2),3));
        dat.alphaCurveFrontal{gi}(:,cond_i,:) = squeeze(mean(mean(gndTF{gi,cond_i}.indv(:,chanID.frontal(log1),freqID.alpha,timeID),2),3));
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

clear Xa
Xa(:,1) = [reshape(dat.betaAvgOccip{1},size(dat.betaAvgOccip{1},1)*size(dat.betaAvgOccip{1},2),1);...
    reshape(dat.betaAvgOccip{2},size(dat.betaAvgOccip{2},1)*size(dat.betaAvgOccip{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.betaOccip, DFs.betaOccip, MSQs.betaOccip, Fs.betaOccip, Ps.betaOccip]=mixed_between_within_anova(Xa);

clear Xa
Xa(:,1) = [reshape(dat.alphaAvgOccip{1},size(dat.alphaAvgOccip{1},1)*size(dat.alphaAvgOccip{1},2),1);...
    reshape(dat.alphaAvgOccip{2},size(dat.alphaAvgOccip{2},1)*size(dat.alphaAvgOccip{2},2),1)];
Xa(:,2) = [ones(sum(subs.group==1)*3,1);ones(sum(subs.group==2)*3,1)*2];
Xa(:,3) = [ones(sum(subs.group==1),1);ones(sum(subs.group==1),1)*2;ones(sum(subs.group==1),1)*3;...
    ones(sum(subs.group==2),1);ones(sum(subs.group==2),1)*2;ones(sum(subs.group==2),1)*3];
Xa(:,4) = [[1:sum(subs.group==1) 1:sum(subs.group==1) 1:sum(subs.group==1)],...
    [1:sum(subs.group==2) 1:sum(subs.group==2) 1:sum(subs.group==2)]+sum(subs.group==1)]';

[SSQs.alpha, DFs.alpha, MSQs.alpha, Fs.alpha, Ps.alpha]=mixed_between_within_anova(Xa);


%% bar plot -beta
addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[100 100 700 1000]);

subplot(3,2,5);hold all;axis square
mn = cellfun(@mean,dat.betaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
plot(xcord(1,:),dat.betaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
errorbar(xcord,mn,se,'k.')

legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%d')
ylabel('PLV(% change)')
title(sprintf('Frontal beta(%.1f~%.1fs)',timeROI(1),timeROI(2)))

if Ps.betaFrontal{3}<.05
    for gi = 1:2
        tmpdata = [dat.betaAvgFrontal{gi} dat.betaAvgFrontal{gi}(:,1)];
        tmpdata = diff(tmpdata,1,2);
        [~,tmp_pval]= ttest(tmpdata);
        sigH = linspace(-0.8,-0.6,3);
        xposi = {[1 2],[2 3],[1 3]};
        for ss = 1:3
            if tmp_pval(ss)<=.05
                plot(xcord(gi,xposi{ss}),[1 1]*sigH(ss),'k','HandleVisibility','off')
                text(mean(xcord(gi,xposi{ss})),sigH(ss),'*','FontSize',18,'HorizontalAlignment','center')
            end
        end
    end
end

subplot(3,2,6);hold all;axis square
mn = cellfun(@mean,dat.betaAvgOccip,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.betaAvgOccip,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.betaAvgOccip{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.betaAvgOccip{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%d')
ylabel('PLV(% change)')
title(sprintf('Posterior beta(%.1f~%.1fs)',timeROI(1),timeROI(2)))

myColors = crameri('batlowW',5);

times = gndTF{1}.time(timeID);
for gi = 1:2
    subplot(3,2,gi);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveFrontal{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveFrontal{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%d')
    ylabel('PLV(% change)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Frontal %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))

    subplot(3,2,gi+2);hold all;
    for cond_i = 1:3
        mn = squeeze(mean(dat.betaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.betaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%d')
    ylabel('PLV(% change)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
end

saveas(gca,fullfile(Dir.figs,['PLVprobe_beta_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preprobe+1,4},'.bmp']))
%%
myColors = crameri('grayC',3);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

figure('Position',[100 100 700 1000]);

subplot(3,2,5);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgOccip,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgOccip,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
plot(xcord(1,:),dat.alphaAvgOccip{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgOccip{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%d')
ylabel('PLV(% change)')
title(sprintf('Posterior alpha(%.1f~%.1fs)',timeROI(1),timeROI(2)))

subplot(3,2,6);hold all;axis square
mn = cellfun(@mean,dat.alphaAvgFrontal,'UniformOutput',false);
mn = vertcat(mn{:});
tmp_std = cellfun(@std,dat.alphaAvgFrontal,'UniformOutput',false);
se = vertcat(tmp_std{:})./[sqrt(sum(subs.group==1)) sqrt(sum(subs.group==2))]';

plot(xcord(1,:),dat.alphaAvgFrontal{1},'Color',[0.8 0.8 0.8],'HandleVisibility','off');
plot(xcord(2,:),dat.alphaAvgFrontal{2},'Color',[0.8 0.8 0.8],'HandleVisibility','off');

hb = bar(mn,'FaceAlpha',0.8);
hb(1).FaceColor = myColors(1,:);
hb(2).FaceColor = myColors(2,:);
hb(3).FaceColor = myColors(3,:);
xcord = vertcat(hb(:).XEndPoints)';
errorbar(xcord,mn,se,'k.')
legend(condStr)
set(gca,'xtick',[1 2],'XTickLabel',groupStr)
ytickformat('%d')
ylabel('PLV(% change)')

title(sprintf('Frontal alpha(%.1f~%.1fs)',timeROI(1),timeROI(2)))
% plot significance
% plot(xcord(1,[2 3]),[1.2 1.2],'k','HandleVisibility','off')
% text(mean(xcord(1,[2 3])),1.2,'***','FontSize',18,'HorizontalAlignment','center')
% plot(xcord(1,[1 3]),[1 1],'k','HandleVisibility','off')
% text(mean(xcord(1,[1 3])),1,'***','FontSize',18,'HorizontalAlignment','center')

myColors = crameri('batlowW',5);

times = gndTF{1}.time(timeID);
for gi = 1:2
    subplot(3,2,gi);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveOccip{gi}(:,cond_i,:)));
        se = squeeze(std(dat.alphaCurveOccip{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%d')
    ylabel('PLV(% change)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Occip %d~%dHz',freq.alphaFreq(1),freq.alphaFreq(2)))

    subplot(3,2,gi+2);hold all;axis square;box on
    for cond_i = 1:3
        mn = squeeze(mean(dat.alphaCurveFrontal{gi}(:,cond_i,:)));
        se = squeeze(std(dat.alphaCurveFrontal{gi}(:,cond_i,:),0,1)./sqrt(sum(subs.group==gi)));
        shadedErrorBar(times,mn,se,{'color',myColors(cond_i,:)})
    end

    plot(get(gca,'XLim'),[0 0],'k--','HandleVisibility','off')
    plot([0 0],get(gca,'YLim'),'k','HandleVisibility','off')
    legend(condStr)
    ytickformat('%d')
    ylabel('PLV(% change)');
    xlabel('Time(0s=maintenance)')
    title(groupStr{gi},sprintf('Frontal %d~%dHz',freq.alphaFreq(1),freq.alphaFreq(2)))
end

saveas(gca,fullfile(Dir.figs,['PLVprobe_alpha_',num2str(freq.alphaFreq(1)),'~',num2str(freq.alphaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preprobe+1,4},'.bmp']))

%% correlation
load beha.mat
wanted = ismember(beha.name,subs.name);
beha = beha(wanted,:);
dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));

% dat.beha{1} = table2array(beha(beha.group==1,{'s1_d','s2_d','s3_d'}));
% dat.beha{2} = table2array(beha(beha.group==2,{'s1_d','s2_d','s3_d'}));


addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',4);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

clear g
for gi = 1:2
    [r,pval] = corr(dat.betaAvgFrontal{gi},dat.beha{gi},'type','Spearman');

    for cond = 1:3
        Y = dat.beha{gi}(:,cond);
        X = dat.betaAvgFrontal{gi}(:,cond);

        g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,cond).geom_point(); %Draw a scatter plot
        g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
        g(gi,cond).set_names('x',sprintf('Beta PLV(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2)),'y','Beha'); %Set the title of the axis
        g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
        g(gi,cond).set_title(sprintf('%s SS#%d\nRho = %.3f, p = %.3f',groupStr{gi},cond,r(cond,cond),pval(cond,cond)));
        g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    end

    cond = 4;
    Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
    X = reshape(dat.betaAvgFrontal{gi},size(dat.betaAvgFrontal{gi},1)*size(dat.betaAvgFrontal{gi},2),1);
    [r,pval] = corr(Y,X,'type','Spearman');

    g(gi,cond) = gramm('x',X,'y',Y); %Specify the values of the horizontal axis x and the vertical axis y, and create a gramm drawing object
    g(gi,cond).geom_point(); %Draw a scatter plot
    g(gi,cond).stat_glm(); %Draw lines and confidence intervals based on scatter plots
    g(gi,cond).set_names('x',sprintf('Beta PLV(%d~%dHz)',freq.betaFreq(1),freq.betaFreq(2)),'y','Beha'); %Set the title of the axis
    g(gi,cond).set_color_options('map',myColors(cond,:)); %Set the color of the point
    g(gi,cond).set_title(sprintf('%s collapse\nRho = %.3f, p = %.3f',groupStr{gi},r,pval));
    g(gi,cond).set_text_options('base_size' ,10,'title_scaling' ,1.1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times


end
figure('Position',[100 100 1100 550]);
g.draw();
saveas(gcf,fullfile(Dir.figs,['PLVprobe_Corr_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI(1)),'~',num2str(timeROI(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preprobe+1,4} '.bmp']))
