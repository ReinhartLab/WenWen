clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0,:) = [];
% subs([11 18 22],:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;

addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))
IsOcci = 0;
IsBL2preDelay = 0;

txtCell = {'','','';'_occi','_bl2preDelay','_corrTrials'};

% permFile = ['Perm_SSdelay',txtCell{IsOcci+1,1},txtCell{IsBL2preDelay+1,2},'.mat'];

%%
clear fem tmp_data
for sub_i = 1:subN
    subname = subs.name{sub_i};
    for IsCorretTrials = [0 1]
        load(fullfile(Dir.results,[subname,'_ssDelay',txtCell{IsOcci+1,1},txtCell{IsBL2preDelay+1,2},txtCell{IsCorretTrials+1,3},'_mal.mat']))
        tmp_data(sub_i,IsCorretTrials+1,:) = dcd.ACC;
    end
end

tmp_data = smoothdata(tmp_data,3,'gaussian',10);%subN*cond*time

for g = 1:2
    tmpID = subs.group==g;
    fem.sl.data{g} = tmp_data(tmpID,:,:);%
    fem.sl.mn{g} = squeeze(mean(fem.sl.data{g})); % cond*time
end

timeRange = dcd.Dtime;
times = dcd.times./1000;
nSamps = dcd.nSamps;
%% bootstrap for SEM across subjects
% if exist(permFile,'file')==0
bIter = 1000;

for g = 1:2
    subNg(g) = sum(subs.group==g);
    boot.IDX = nan(bIter,subNg(g)); %bootstrap sample index
    boot.slope = nan(bIter,subNg(g),size(fem.sl.data{g},2),nSamps);%bootstrapped CTFs
    boot.M = nan(bIter,size(fem.sl.data{g},2),nSamps);% the mean osbserved slope (across subs)

    fprintf('Bootstrap replication...\n')
    for b = 1:bIter
        [boot.slope(b,:,:,:), boot.IDX(b,:)] = datasample(fem.sl.data{g},subNg(g),1); % sample nSubs many observations for the subs dimensions (with replacement)
        boot.M(b,:,:) = mean(boot.slope(b,:,:,:));
    end
    boot.SEM = squeeze(std(boot.M));
    boot = rmfield( boot , {'slope'});
    fem.boot{g} = boot;clear boot

    for cond = 1:2
        dat = squeeze(fem.sl.data{g}(:,cond,:))'-1/3;%time*sub
        diffstat = 'diff';%'t';% or 'diff'
        nperm = 10000;

        [datobs, datrnd] = cluster_test_helper(dat, nperm, diffstat);
        tail = 0;alpha = .05;
        clusteralpha = .05;
        clusterstat = 'sum'; %'size'

        [perm(g,cond).pval, perm(g,cond).h, perm(g,cond).clusterinfo] = cluster_test(datobs, datrnd, tail, alpha,...
            clusteralpha, clusterstat);
    end
end
% save(permFile,'fem','perm','-v7.3')

% else
%     load(permFile)
% end
%%
condStr = {'All trials','Correct trials'};
groupStr = {'Young','Old'};

myFigBasic;

myColors = crameri('bam',2);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

lineW = 1.2;

figure('name',['SSdelay',txtCell{IsOcci+1,1},txtCell{IsBL2preDelay+1,2}] ,'position',[100 100 900 300]);

sigValue = linspace(0.31,0.32,2);

for cond = 1:2
    subplot(1,2,cond);
    hold all;axis square
    for g = 1:2
        shadedErrorBar(times,fem.sl.mn{g}(cond,:),fem.boot{g}.SEM(cond,:),{'color',myColors(g,:),'LineWidth',lineW});

        sig =  nan(size(perm(g,cond).h));
        sig(perm(g,cond).h) = sigValue(g);
        plot(times,sig,'color', myColors(g,:),'LineWidth',2,'HandleVisibility','off')
    end

    % sig =  nan(size(perm(1).dif.h));
    % sig(perm(1).dif.h) =sigValue(3);
    % plot(times,sig,'k','LineWidth',2,'HandleVisibility','off')
    ytickformat('%,.2f')

    xlabel({'\bfTime (s)'});ylabel('\bfAccuracy');
    legend(groupStr,'Location','bestoutside');
    set(gca,'xtick',0:1:3,'xlim',dcd.Dtime,'ytick',0.3:0.05:0.65,'ylim',[0.3 0.5]);
    plot(get(gca,'xlim'),[1/3 1/3],'k','HandleVisibility','off')
    plot([0 0],get(gca,'ylim'),'k','HandleVisibility','off')
    plot([3 3],get(gca,'ylim'),'k','HandleVisibility','off')
    plot([3.2 3.2],get(gca,'ylim'),'k:','HandleVisibility','off')
    title(condStr{cond});
    hold off
end
saveas(gcf,fullfile(Dir.figs,['DecodeSS' txtCell{IsOcci+1,1},txtCell{IsBL2preDelay+1,2} '.bmp']))