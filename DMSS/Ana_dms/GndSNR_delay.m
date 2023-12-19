clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
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

% timeROI.all = [-5.3 4.5];% in s
% SStime = {[-1.6 timeROI.all(2)],[-2.8 timeROI.all(2)],[-5.2 timeROI.all(2)]};% epoch length for different SS

timeROI.delay = [0 3];% in s

%%
clear subsAll
for sub_i = 1:subN
    subname = subs.name{sub_i};
    matName = fullfile(Dir.results,[subname,'_SNR_delay_',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);
        
        subsAll.broad(sub_i,:,:,:) = SNR.broad;

        subsAll.freq(sub_i,:,:,:,:) = SNR.freqdata;
    end
end

%%
load chanLocs_dms.mat
chanloc(contains({chanloc.labels},'EOG')) = [];
[~,chanID] = ismember(chanROI,{chanloc.labels});

freqID.beta = dsearchn(SNR.freq',freq.betaFreq');
freqID.beta = freqID.beta(1):freqID.beta(2);

timeID.delay = dsearchn(SNR.times',timeROI.delay');
timeID.delay = timeID.delay(1):timeID.delay(2);

subsAll.avg.broad = squeeze(mean(mean(subsAll.broad(:,:,chanID,timeID.delay),3),4));%subs*ss

timeID.toi = dsearchn(SNR.ftoi',timeROI.delay');
timeID.toi = timeID.toi(1):timeID.toi(2);

subsAll.avg.freq = squeeze(mean(mean(subsAll.freq(:,:,chanID,:,timeID.toi),3),5));%subs*ss*freq
subsAll.avg.beta = squeeze(mean(mean(mean(subsAll.freq(:,:,chanID,freqID.beta,timeID.toi),3),4,'omitnan'),5,'omitnan'));%subs*ss

SNR = repmat(subs(:,["name","group","groupStr"]),3,1);
SNR.ss = sort(repmat(1:3,1,subN))';
SNR.ss = categorical(SNR.ss);
SNR.broad = reshape(subsAll.avg.broad,subN*3,1);
SNR.freq = reshape(subsAll.avg.beta,subN*3,1);

%% anova
Xa = SNR(:,["broad","group","ss"]);
Xa.ss = double(Xa.ss);
Xa = table2array(Xa);
Xa(:,4) = repmat(1:subN,1,3)';
[SSQs.broad, DFs.broad, MSQs.broad, Fs.broad, Ps.broad]=mixed_between_within_anova(Xa);

Xa = SNR(:,["freq","group","ss"]);
Xa.ss = double(Xa.ss);
Xa = table2array(Xa);
Xa(:,4) = repmat(1:subN,1,3)';
[SSQs.freq, DFs.freq, MSQs.freq, Fs.freq, Ps.freq]=mixed_between_within_anova(Xa);


%%
figure('Position',[100 100 500 300]);
subplot(1,2,1);hold all;axis square
boxchart(SNR.ss,SNR.broad,'GroupByColor',SNR.groupStr)
legend('Location','southoutside')
xlabel('set sizes')
ylabel('SNR=mean/sd','Interpreter','none')
text(0.1,0.9,sprintf('Betw: F = %.3f, p = %.3f\nWith: F = %.3f, p = %.3f\nInter: F = %.3f, p = %.3f',Fs.broad{1},Ps.broad{1},Fs.broad{3},Ps.broad{3},Fs.broad{4},Ps.broad{4}),'sc')
title('Broad band')

subplot(1,2,2);hold all;axis square
boxchart(SNR.ss,SNR.freq,'GroupByColor',SNR.groupStr)
legend('Location','southoutside')
xlabel('set sizes')
ylabel('SNR=mean/sd','Interpreter','none')
text(0.1,0.9,sprintf('Betw: F = %.3f, p = %.3f\nWith: F = %.3f, p = %.3f\nInter: F = %.3f, p = %.3f',Fs.freq{1},Ps.freq{1},Fs.freq{3},Ps.freq{3},Fs.freq{4},Ps.freq{4}),'sc')
title(sprintf('%d~%dHz',freq.betaFreq(1),freq.betaFreq(2)))
saveas(gcf,fullfile(Dir.figs,['SNR_delay_',[chanROI{:}],num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz_',txtCell{IsLap+1,1},num2str(timeROI.delay(1)),'~',num2str(timeROI.delay(2)),'s',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
