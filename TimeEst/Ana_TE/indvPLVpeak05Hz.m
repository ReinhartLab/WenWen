clear
load('subs.mat');
txtCell = {'','','';'_lap','_dephase','_bl'};
IsLap = 0;
IsdePhase = 1;
IsBL = 1;% whether to apply baseline

condStr = {'Correct','Incorrect','Incor-Corr'};
PeakPara.time = [0.1 0.5]; % in s
PeakPara.freq = [3 12]; % in Hz
% PeakPara.freq = [3 10]; % in Hz

%%
for sn = 1:height(subs)
    subname = subs.name{sn};
    if subs.excluded(sn)==1
        continue
    end
    Indv.subname{sn,1} = subname;

    outFile = fullfile(Dir.results,[subname,'_con0.5Hz',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},'.mat']);
    if isfile(outFile)
        load(outFile)
    end

    if IsBL % baseline correction
        for cond_i = 1:2
            tfDat{cond_i}.plvFT.plvspctrm = bsxfun(@minus,tfDat{cond_i}.plvFT.plvspctrm,mean(tfDat{cond_i}.plvBS.plvspctrm,3));
        end
    end

    for cond_i = 1:2
        tfAll.subs(sn,cond_i,:,:,:) = tfDat{cond_i}.plvFT.plvspctrm;
    end

    cond_i = 3;%incor-corr
    tfDat{cond_i}.plvFT = tfDat{1}.plvFT;
    tfDat{cond_i}.plvFT.plvspctrm = tfDat{2}.plvFT.plvspctrm - tfDat{1}.plvFT.plvspctrm;
    tfAll.subs(sn,cond_i,:,:,:) = tfDat{cond_i}.plvFT.plvspctrm;
    %% time freq map
    myFigBasic
    w = [1 1 1];
    b = [0,0,1];
    k = [0,0,0];
    r = [1,0,0];
    bkr = createcolormap(w,b,k,r,w); % 256x3 array

    cfg = [];
    cfg.parameter = 'plvspctrm';
    cfg.colormap = bkr;
    cfg.refchannel = {'FCz'};
    cfg.channel = {'F6'};
    cfg.figure = 'gca';

    figure('Position',[100 100 1000 600]);

    for cond_i = 1:3
        subplot(2,3,cond_i);hold all;axis square
        ft_singleplotTFR(cfg, tfDat{cond_i}.plvFT);

        title(condStr{cond_i},sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]));colorbar
        ylabel('\bfFrequency(Hz)');
        xlabel('\bfTime(s)')
    end

    tmpfreq = tfDat{1}.plvFT.freq<=PeakPara.freq(2) & tfDat{1}.plvFT.freq>=PeakPara.freq(1);
    tmptoi = tfDat{1}.plvFT.time<=PeakPara.time(2) & tfDat{1}.plvFT.time>=PeakPara.time(1);
    tmpChancomb = ismember(tfDat{1}.plvFT.labelcmb(:,1),cfg.refchannel) & ismember(tfDat{1}.plvFT.labelcmb(:,2),cfg.channel);

    for cond_i = 1:3
        freqCurve(cond_i,:) = squeeze(mean(tfDat{cond_i}.plvFT.plvspctrm(tmpChancomb,:,tmptoi),3));% average time
    end

    tmpdat = freqCurve(3,:);%use difference for peak detection
    [pks,locs] = findpeaks(tmpdat);
    [pks2,locs2] = max(tmpdat(tmpfreq));% in case no peaks detected
    locs2 = locs2 + dsearchn(tfDat{1}.plvFT.freq',PeakPara.freq(1))-1;
    pks = unique([pks pks2],'stable');% combine both methods
    locs = unique([locs locs2],'stable');

    tmp = tfDat{1}.plvFT.freq(locs)<=PeakPara.freq(2) & tfDat{1}.plvFT.freq(locs)>=PeakPara.freq(1);
    Indv.peakF{sn,1} = round(tfDat{1}.plvFT.freq(locs(tmp)),1);

    subplot(2,1,2);hold all;box on
    plot(tfDat{1}.plvFT.freq,freqCurve)
    plot(tfDat{1}.plvFT.freq(locs),pks,'ro')
    for li = 1:length(pks)
        text(tfDat{1}.plvFT.freq(locs(li)),pks(li)+0.02,sprintf('%.1fHz',tfDat{1}.plvFT.freq(locs(li))),'HorizontalAlignment','center')
    end
    legend(condStr,'Location','eastoutside')
    xlabel('Frequency(Hz)');
    ylabel('PLV')
    title(sprintf('sn%02d(%s):%.1f~%.1fs',sn,subname,PeakPara.time(1),PeakPara.time(2)));

    saveas(gcf,fullfile(Dir.figs,'IndvPLVpeak0.5Hz',[sprintf('sn%02d_%s%.1f~%.1fs',sn,subname,PeakPara.time(1),PeakPara.time(2)),sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.png']))
end
save(fullfile(Dir.results,['indvPeakF0.5Hz.mat']),"Indv")
%% peakF distribution
load(fullfile(Dir.results,['indvPeakF0.5Hz.mat']))
Indv = struct2table(Indv);
Indv(subs.excluded==1,:) = [];
Indv.peakF2plot = Indv.peakF;

for i = 1:height(Indv)
    if length(Indv.peakF2plot{i})>1
        Indv.peakF2plot{i} =  min(Indv.peakF2plot{i});
    end
end

myFigBasic;
addpath(genpath('D:\intWM-E\toolbox\RainCloudPlots-master'))
figure('Position',[200 200 160 300]);
hold all;box on
peakF = cell2mat(Indv.peakF2plot);
dodge_shift = 0.12;
[h, u] = raincloud_plot(peakF,'box_on', 1,'color',[0 0 0],...
    'alpha', 0.3,'box_dodge', 1, 'box_dodge_amount', -dodge_shift, 'dot_dodge_amount', -dodge_shift);
xlabel('Peak frequency(Hz)')

set(gca,'xtick',2:4:16, 'YLim', [0 0.25],'ytick',[0 0.25]);
ylabel({'\bfDensity'});
view(90, -90)
set(gca,'TickLength', [0 0]);
text(0.1,0.9,sprintf('Mean = %.1f\nMedian =%.1f',mean(peakF),median(peakF)),'Units','normalized')
saveas(gcf,fullfile(Dir.figs,'IndvPLVpeak0.5Hz',[sprintf('Boxplot%.1f~%.1fs%d~%dHz',PeakPara.time(1),PeakPara.time(2),PeakPara.freq(1),PeakPara.freq(2)),sprintf('%s (ref at %s)',[cfg.channel{:}],[cfg.refchannel{:}]),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL+1,3},'.png']))
