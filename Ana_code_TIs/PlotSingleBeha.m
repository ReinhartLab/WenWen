function PlotSingleBeha(subname)

load subsinfo.mat
load(fullfile(Dir.beha.IBL,subname,[subname,'.mat']));
clear ACC RT
ACC.avgTask = nan(3,6,params.trlNperTaskperRun);% three stim phases * 6 blocks* repeat
RT.avgTask = ACC.avgTask;
for p = 1:2
    for bi = 1:6
        tmp_table = subsMat(p).table(bi).all;
        for ri = 1:params.trlNperTaskperRun
            ACC.avgTask(p,bi,ri) = sum(tmp_table.acc(tmp_table.repi == ri)==1)/3;% averaged acc across three tasks
            if ACC.avgTask(p,bi,ri)~=0
                RT.avgTask(p,bi,ri) = mean(tmp_table.RT(tmp_table.repi == ri & tmp_table.acc == 1),'omitnan');
            end
        end

        for t = 1:3
            tmpIdx = tmp_table.curTask==t;
            ACC.eachTask(p,bi,t) = sum(tmp_table.acc(tmpIdx)==1)/params.trlNperTaskperRun;
            RT.eachTask(p,bi,t) = mean(tmp_table.RT(tmp_table.acc==1 & tmpIdx));
        end
    end
end

%%  plot
myColors = jet(6);
taskStr = params.taskStr;
titleStr = {'Pre','Stim'};
blockStr = {'Block1','Block2','Block3','Block4','Block5','Block6'};
figure('Name',subname,'Position',[0 100 1500 500]);
for p = 1:2
    subplot(1,5,p);
    hold all;axis square
    for bi = 1:6
        yyaxis left
        plot(squeeze(RT.avgTask(p,bi,:)),'-','Color',myColors(bi,:),'LineStyle','-')

        yyaxis right
        plot(squeeze(ACC.avgTask(p,bi,:)),'o--','Color',myColors(bi,:),'HandleVisibility','off')
        set(gca,'ylim',[0 1])
    end
    yyaxis left
    plot(squeeze(mean(RT.avgTask(p,:,:),2,'omitnan')),'k-','LineWidth',3,'HandleVisibility','off')
    ylabel('RT.avgTask(s)')
    yyaxis right
    plot(squeeze(mean(ACC.avgTask(p,:,:),2)),'ko--','LineWidth',3,'HandleVisibility','off')
    xlabel('Trial')
    ylabel('Accuracy');

    legend(blockStr,'Location','southoutside')
    title([titleStr{p},num2str(subInfo.stim)])
end

subplot(1,5,3);hold all;axis square
for p = 1:2
    y = squeeze(mean(RT.avgTask(p,:,:),2,'omitnan'));% array of tmp_table
    x = 1:6;
    slope(p,:) = polyfit(x,y,1);% slope and intercept
end
plot(slope(:,1),'k-o')
title('slope')
set(gca,'xticklabel',titleStr,'XLim',[1 2])

subplot(1,5,4);hold all;axis square
bar(squeeze(mean(ACC.eachTask,2,'omitnan')));% si*bi*ti
legend(taskStr,'Location','southoutside')
title('Accuracy')
set(gca,'xticklabel',titleStr,'Xtick',1:2,'ytick',[0.4:0.2:1],'YLim',[0.4 1])

subplot(1,5,5);hold all;axis square
bar(squeeze(mean(RT.eachTask,2,'omitnan')));% si*bi*ti
legend(taskStr,'Location','southoutside')
title('RT')
set(gca,'xticklabel',titleStr,'Xtick',1:2)
sgtitle(subname)

saveas(gca,fullfile(Dir.figs,'BehaIndv',[subname '.png']))
