clear;rng('shuffle');

for IsBL2preDelay = [0 1]
    txtCell = {'','';'_bl2preDelay','_corrTrials'};
    load('subs.mat');

    subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
    subN = height(subs);
    %%
    clear tmp_data
    for sub_i = 1:subN
        subname = subs.name{sub_i};
        for IsCorretTrials = [0 1]
            matName =fullfile(Dir.results,[subname,'_ssDelayChan',txtCell{IsBL2preDelay+1,1},txtCell{IsCorretTrials+1,2},'_mal.mat']);

            if isfile(matName)
                load(matName);
                tmp_data(sub_i,IsCorretTrials+1,:,:) = squeeze(dcd.ACC);%2*chan*time
            end

        end
    end
    empty_idx = tmp_data(:,1,1,1)==0;
    tmp_data(empty_idx,:,:,:) = [];
    subs(empty_idx,:) = [];

    tmp_data = smoothdata(tmp_data,4,'gaussian',10);%subN*2*chan*time

    for gi = 1:2
        tmpID = subs.group==gi;
        fem.sl.data{gi} = tmp_data(tmpID,:,:,:);%
        fem.sl.mn{gi} = squeeze(mean(fem.sl.data{gi},'omitnan')); % 2*chan*time
    end

    timeRange = dcd.Dtime;
    times = dcd.times./1000;
    nSamps = dcd.nSamps;

    %%
    condStr = {'All trials','Correct trials'};
    groupStr = {'Young','Old','Y minus O'};
    toi = -1:0.5:3.2;% in s

    [N,edges] = discretize(times,toi,'IncludedEdge','right');

    myFigBasic;

    load chanLocs_dms.mat
    [is, loc] = ismember({dcd.neighbours.label},{chanloc.labels});
    chanloc  = chanloc(loc(is));

    threshP = 0.005;
    threshPdiff = 0.005;
    for cond = 1:2
        figure('name',['SSdelaySearchLight' txtCell{IsBL2preDelay+1,1}],'Position',[10 10 2000 600]);

        for gi = 1:2
            for t = 1:length(toi)-1
                subplot(3,length(toi),(gi-1)*length(toi)+t)

                dat_mn = squeeze(mean(fem.sl.mn{gi}(cond,:,N==t),3,"omitnan"));
                dat_all = squeeze(mean(fem.sl.data{gi}(:,cond,:,N==t),4,'omitnan'));
                [~,p] = ttest(dat_all,1/3,'Tail','right');

                topoplot(dat_mn,chanloc,'emarker2',{find(p<threshP),'o','k',4},'plotrad',0.53,'electrodes','off' );
                caxis([0.32, 0.36])
                h = colorbar;
                %             set(h,'ytick',[0.33:0.03:0.4])
                h.Label.String = 'Classification accuracy';
                h.Label.Rotation = -90;
                h.Label.Position = h.Label.Position + [1 0 0];
                title(sprintf('%s:%.1f~%.1fs', groupStr{gi},toi(:,t),toi(:,t+1)),sprintf('%s, p<%.3f',condStr{cond},threshP));
            end

            t = length(toi); % last column, avg
            subplot(3,length(toi),(gi-1)*length(toi)+t)
            hold all;axis square

            dat_mn = squeeze(mean(fem.sl.mn{gi}(cond,:,:),3));
            dat_all = squeeze(mean(fem.sl.data{gi}(:,cond,:,:),4));
            [~,p] = ttest(dat_all,1/3,'Tail','right');

            topoplot(dat_mn,chanloc,'emarker2',{find(p<threshP),'o','k',4},'plotrad',0.53,'electrodes','off' );
            caxis([0.32, 0.38])
            h = colorbar;
            %             set(h,'ytick',[0.33:0.03:0.4])
            h.Label.String = 'Classification accuracy';
            h.Label.Rotation = -90;
            h.Label.Position = h.Label.Position + [1 0 0];

            title(sprintf('%s:0~3s', groupStr{gi}),sprintf('%s, p<%.3f',condStr{cond},threshP));

        end
        gi = 3;% diff
        for t = 1:length(toi)-1

            subplot(3,length(toi),(gi-1)*length(toi)+t)
            hold all;axis square

            tmp1 = mean(fem.sl.mn{1}(cond,:,N==t),3);
            tmp2 = mean(fem.sl.mn{2}(cond,:,N==t),3);
            dat_mn = squeeze(tmp1-tmp2);
            dat_all1 = squeeze(mean(fem.sl.data{1}(:,cond,:,N==t),4));
            dat_all2 = squeeze(mean(fem.sl.data{2}(:,cond,:,N==t),4));
            [~,p] = ttest2(dat_all1,dat_all2);

            topoplot(dat_mn,chanloc,'emarker2',{find(p<threshPdiff),'o','k',4},'plotrad',0.53,'electrodes','off');
            caxis([-0.04, 0.04])

            h = colorbar;
            %             set(h,'ytick',[-0.02 0.02])
            title(sprintf('%s:%s (%.1fs)', groupStr{gi},condStr{cond},mean(toi(:,t))),sprintf('p<%.2f',threshPdiff));
        end

        t = length(toi); % last column, avg

        subplot(3,length(toi),(gi-1)*length(toi)+t)
        hold all;axis square

        tmp1 = mean(fem.sl.mn{1}(cond,:,:),3);
        tmp2 = mean(fem.sl.mn{2}(cond,:,:),3);
        dat_mn = squeeze(tmp1-tmp2);
        dat_all1 = squeeze(mean(fem.sl.data{1}(:,cond,:,:),4));
        dat_all2 = squeeze(mean(fem.sl.data{2}(:,cond,:,:),4));
        [~,p] = ttest2(dat_all1,dat_all2);

        topoplot(dat_mn,chanloc,'emarker2',{find(p<threshPdiff),'o','k',4},'plotrad',0.53,'electrodes','off');
        caxis([-0.04, 0.04])

        h = colorbar;
        %         set(h,'ytick',[-0.02 0.02])
        title(sprintf('%s:%s avg', groupStr{gi},condStr{cond}),sprintf('p<%.3f',threshPdiff));

        saveas(gcf,fullfile(Dir.figs,['DecodeSSdelaySearchLight_',condStr{cond},txtCell{IsBL2preDelay+1,1},'.bmp']))
    end
end