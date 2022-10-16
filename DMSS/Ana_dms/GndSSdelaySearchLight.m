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
groupStr = {'Young','Old','O minus Y'};

toi = 0:0.5:3;% in s
tw = 50;%in ms
toi = toi+[-1;1]*tw/1000 ;
clear topoIdx
topoIdx(:,1) = dsearchn(times',toi(1,:)');
topoIdx(:,2) = dsearchn(times',toi(2,:)');

myFigBasic;

load chanLocs_dms.mat
[is, loc] = ismember({dcd.neighbours.label},{chanloc.labels});
chanloc  = chanloc(loc(is));

for cond = 1:2
    ftopo= figure('name',['SSdelaySearchLight' txtCell{IsBL2preDelay+1,1}],'Position',[10 10 2000 500]);

    for g = 1:2
        for t = 1:length(toi)
            subplot(3,length(toi),(g-1)*length(toi)+t)
            hold all;axis square

            dat = squeeze(mean(fem.sl.mn{g}(cond,:,topoIdx(t,1):topoIdx(t,2)),3));

            topoplot(dat,chanloc);
            caxis([0.32, 0.38])
            h = colorbar;
            set(h,'ytick',[0.32 0.38])
            %             h.Label.String = 'acc';
            %             h.Limits = [0.3, 0.4];
            title(sprintf('%s:%s (%.2fs)', groupStr{g},condStr{cond},mean(toi(:,t))));
        end
    end
    g = 3;% old minus young
    for t = 1:length(toi)

        subplot(3,length(toi),(g-1)*length(toi)+t)
        hold all;axis square

        tmp1 = mean(fem.sl.mn{1}(cond,:,topoIdx(t,1):topoIdx(t,2)),3);
        tmp2 = mean(fem.sl.mn{2}(cond,:,topoIdx(t,1):topoIdx(t,2)),3);
        dat = squeeze(tmp2-tmp1);

        topoplot(dat,chanloc);
        caxis([-0.02, 0.02])

        h = colorbar;
        set(h,'ytick',[-0.02 0.02])
        title(sprintf('%s:%s (%.2fs)', groupStr{g},condStr{cond},mean(toi(:,t))));
    end

    saveas(gcf,fullfile(Dir.figs,['SSdelaySearchLight',condStr{cond},txtCell{IsBL2preDelay+1,1},'.bmp']))
end