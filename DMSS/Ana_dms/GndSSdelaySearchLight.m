clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\Toolbox\fieldtrip-20211209')
ft_defaults;

IsBL2preDelay = 0;
txtCell = {'','';'_bl2preDelay','_corrTrials'};

%%
clear fem tmp_data
for sub_i = 1:subN
    subname = subs.name{sub_i};
    for IsCorretTrials = [0 1]
        load(fullfile(Dir.results,[subname,'_ssDelayChan',txtCell{IsBL2preDelay+1,1},txtCell{IsCorretTrials+1,2},'_mal.mat']))

        tmp_data(sub_i,IsCorretTrials+1,:,:) = squeeze(dcd.ACC);%2*chan*time
    end
end

tmp_data = smoothdata(tmp_data,4,'gaussian',10);%subN*2*chan*time

for g = 1:2
    tmpID = subs.group==g;
    fem.sl.data{g} = tmp_data(tmpID,:,:,:);%
    fem.sl.mn{g} = squeeze(mean(fem.sl.data{g})); % 2*chan*time
end

timeRange = dcd.Dtime;
times = dcd.times./1000;
nSamps = dcd.nSamps;

%%
condStr = {'All','Correct'};
groupStr = {'Young','Old','Y minus O'};
toi = 0:0.4:3.2;% in s

N = discretize(times,toi,'IncludedEdge','right');

myFigBasic;

load chanLocs_dms.mat
[is, loc] = ismember({dcd.neighbours.label},{chanloc.labels});
chanloc  = chanloc(loc(is));

threshP = 0.001;
threshPdiff = 0.05;
for cond = 1:2
    ftopo= figure('name',['SSdelaySearchLight' txtCell{IsBL2preDelay+1,1}],'Position',[10 10 2000 600]);

    for g = 1:2
        for t = 1:length(toi)-1
            subplot(3,length(toi),(g-1)*length(toi)+t)
            hold all;axis square

            dat_mn = squeeze(mean(fem.sl.mn{g}(cond,:,N==t),3));
            dat_all = squeeze(mean(fem.sl.data{g}(:,cond,:,N==t),4));
            [~,p] = ttest(dat_all,1/3,'Tail','right');

            topoplot(dat_mn,chanloc,'emarker2',{find(p<threshP),'o','k',4},'plotrad',0.53,'electrodes','off' );
            caxis([0.31, 0.38])
            h = colorbar;
            set(h,'ytick',[0.33:0.03:0.4])
            %             h.Label.String = 'acc';
            %             h.Limits = [0.3, 0.4];
            title(sprintf('%s:%s (%.1fs)', groupStr{g},condStr{cond},mean(toi(:,t))),sprintf('p<%.3f',threshP));
        end
    end
    g = 3;% old minus young
    for t = 1:length(toi)-1

        subplot(3,length(toi),(g-1)*length(toi)+t)
        hold all;axis square

        tmp1 = mean(fem.sl.mn{1}(cond,:,N==t),3);
        tmp2 = mean(fem.sl.mn{2}(cond,:,N==t),3);
        dat_mn = squeeze(tmp1-tmp2);
        dat_all1 = squeeze(mean(fem.sl.data{1}(:,cond,:,N==t),4));
        dat_all2 = squeeze(mean(fem.sl.data{2}(:,cond,:,N==t),4));
        [~,p] = ttest2(dat_all1,dat_all2);

        topoplot(dat_mn,chanloc,'emarker2',{find(p<threshPdiff),'o','k',4},'plotrad',0.53,'electrodes','off');
        caxis([-0.02, 0.02])

        h = colorbar;
        set(h,'ytick',[-0.02 0.02])
        title(sprintf('%s:%s (%.2fs)', groupStr{g},condStr{cond},mean(toi(:,t))),sprintf('p<%.2f',threshPdiff));
    end

    saveas(gcf,fullfile(Dir.figs,['SSdelaySearchLight',condStr{cond},txtCell{IsBL2preDelay+1,1},'.bmp']))
end