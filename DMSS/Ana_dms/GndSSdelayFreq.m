clear;rng('shuffle');

addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsOcci = 0;
IsdePhase = 1;
for IsBL2preDelay = [0 1]

    txtCell = {'','','','';'_occi','_dephase','_bl2preDelay','_corrTrials'};
    % permFile = ['Perm_SSdelayFreq',txtCell{IsOcci+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL2preDelay+1,3},'.mat'];
    load('subs.mat');

    subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
    subN = height(subs);
    %%
    clear fem tmp_data
    for sub_i = 1:subN
        subname = subs.name{sub_i};
        for IsCorretTrials = [0 1]
            matName = fullfile(Dir.results,[subname,'_ssDelayFreq',txtCell{IsOcci+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL2preDelay+1,3},txtCell{IsCorretTrials+1,4},'_mal.mat']);

            if isfile(matName)
                load(matName);
                tmp = vertcat(dcd.ACC{:});

                tmp_data(sub_i,IsCorretTrials+1,:,:) = tmp;
            end
        end
    end

    empty_idx = tmp_data(:,1,1,1)==0;
    tmp_data(empty_idx,:,:,:) = [];
    subs(empty_idx,:) = [];
    % tmp_data = smoothdata(tmp_data,4,'gaussian',10);%subN*cond*freq*time

    for gi = 1:2
        tmpID = subs.group==gi;
        fem.sl.data{gi} = tmp_data(tmpID,:,:,:);%subN*2*freq*time
        fem.sl.mn{gi} = squeeze(mean(fem.sl.data{gi})); % cond*freq*time
    end

    timeRange = dcd.Dtime;
    times = dcd.times./1000;
    nSamps = dcd.nSamps;
    freqs = mean(dcd.freqs);

    %% bootstrap for SEM across subjects
    % if exist(permFile,'file')==0
    bIter = 1000;

    for gi = 1:2
        %     subNg(g) = sum(subs.group==g);
        %     boot.IDX = nan(bIter,subNg(g)); %bootstrap sample index
        %     boot.slope = nan(bIter,subNg(g),size(fem.sl.data{g},2),size(fem.sl.data{g},3),nSamps);%bootstrapped CTFs
        %     boot.M = nan(bIter,size(fem.sl.data{g},2),size(fem.sl.data{g},3),nSamps);% the mean osbserved slope (across subs)
        %
        %     fprintf('Bootstrap replication...\n')
        %     for b = 1:bIter
        %         [boot.slope(b,:,:,:,:), boot.IDX(b,:)] = datasample(fem.sl.data{g},subNg(g),1); % sample nSubs many observations for the subs dimensions (with replacement)
        %         boot.M(b,:,:,:) = mean(boot.slope(b,:,:,:,:));
        %     end
        %     boot.SEM = squeeze(std(boot.M));
        %     boot = rmfield( boot , {'slope'});
        %     fem.boot{g} = boot;clear boot

        %     for cond = 1:2
        %         dat = squeeze(fem.sl.data{g}(:,cond,:))'-1/3;%time*sub
        %         diffstat = 'diff';%'t';% or 'diff'
        %         nperm = 10000;
        %
        %         [datobs, datrnd] = cluster_test_helper(dat, nperm, diffstat);
        %         tail = 0;alpha = .05;
        %         clusteralpha = .05;
        %         clusterstat = 'sum'; %'size'
        %
        %         [perm(g,cond).pval, perm(g,cond).h, perm(g,cond).clusterinfo] = cluster_test(datobs, datrnd, tail, alpha,...
        %             clusteralpha, clusterstat);
        %     end
    end
    % save(permFile,'fem','perm','-v7.3')

    % else
    %     load(permFile)
    % end
    %%
    condStr = {'All trials','Correct trials'};
    groupStr = {'Young','Old'};
    myFigBasic;

    lineW = 1.2;

    figure('name',['Perm_SSdelay',txtCell{IsOcci+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL2preDelay+1,3}],'position',[100 100 1000 400]);

    sigValue = linspace(0.31,0.32,2);

    for cond = 1:2%
        for gi = 1:2
            subplot(2,3,(cond-1)*3+gi);
            hold all;axis square
            dat = squeeze(fem.sl.mn{gi}(cond,:,:));
            imagesc(times,freqs,dat);
            colormap jet

            hc = colorbar();
            hc.Label.String = 'acc';
            hc.Limits = [0.32 0.36];

            plot([0 0],get(gca,'ylim'),'k--','HandleVisibility','off')
            plot([3 3],get(gca,'ylim'),'k--','HandleVisibility','off')
            xlabel({'\bfTime (s)'});ylabel('\bfFrequency(Hz)');
            set(gca,'xtick',[0 3],'ytick',[1:4:20],'YTick',[0:10:40]);
            title([condStr{cond} '-' groupStr{gi}]);
        end
    end

    for cond = 1:2
        subplot(2,3,(cond-1)*3+3);

        hold all;axis square
        dat = squeeze(fem.sl.mn{2}(cond,:,:)-fem.sl.mn{1}(cond,:,:));
        imagesc(times,freqs,dat);

        caxis([-0.06, 0.06])

        hc = colorbar();
        hc.Label.String = 'acc';
        hc.Limits = [-0.06, 0.06];

        plot([-1 -1],get(gca,'ylim'),'k--','HandleVisibility','off')
        plot([0 0],get(gca,'ylim'),'k--','HandleVisibility','off')

        plot([3 3],get(gca,'ylim'),'k--','HandleVisibility','off')
        xlabel({'\bfTime (s)'});ylabel('\bfFrequency(Hz)');
        set(gca,'xtick',[0 3],'ytick',[1:4:20],'YTick',[0:10:40]);
        title([condStr{cond} ': old minus young']);

    end
    saveas(gcf,fullfile(Dir.figs,['DecodeFreq',txtCell{IsOcci+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL2preDelay+1,3} ,'.bmp']))
end