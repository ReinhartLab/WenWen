clear;
load subs.mat
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

subs(subs.rawEEG==0 | subs.exclude ==1,:) = [];

subN = height(subs);

%%
outA = 3;
RTmax = 5;% RT limit
clear subsBeha

for sub_i = 1:subN
    csvFile = fullfile(Dir.beha,subs.csvFile{sub_i});

    M = readtable(csvFile);    
    M = M(:,["block_num","ss_num","button_resp_rt","button_resp_corr"]);

    M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});
      
    M(isnan(M.ss_num),:) = [];
    if sum(isnan(M.button_resp_rt)) >0
        M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
    end

%     ssValue = [1 2 4];
%     for ss = 1:3
% 
%         nTarget = sum(ismember(M.type,'d') & M.ss_num == ssValue(ss));
%         nDistract = sum(ismember(M.type,'s') & M.ss_num == ssValue(ss));
% 
%         hit = sum(M.button_resp_corr==1 & ismember(M.type,'d') & M.ss_num == ssValue(ss))/nTarget;
%         FA = sum(M.button_resp_corr==0 & ismember(M.type,'s') & M.ss_num == ssValue(ss))/nDistract;
% 
%         subsBeha(sub_i,ss,3) =  dprimeCalculator(hit,FA,nTarget,nDistract);
%     end

    tmp_sum = groupsummary(M,"ss_num","mean",["button_resp_corr"]);
    subsBeha(sub_i,:,1) = table2array(tmp_sum(:,["mean_button_resp_corr"]));%sub*set*(acc)

    tmp_rt = datastats(M.button_resp_rt(M.button_resp_corr==1 & M.button_resp_rt<RTmax));
    downLim = tmp_rt.mean-outA*tmp_rt.std;
    upLim = tmp_rt.mean+outA*tmp_rt.std;

    tmp_idx = M.button_resp_rt<upLim &  M.button_resp_rt>downLim &  M.button_resp_corr==1 & M.button_resp_rt<RTmax;

    tmp_sum = groupsummary(M(tmp_idx,:),"ss_num","mean",["button_resp_rt"]);
    subsBeha(sub_i,:,2) = table2array(tmp_sum(:,["mean_button_resp_rt"]));%sub*set*(rt)
end

subs.s1_acc = subsBeha(:,1,1);
subs.s2_acc = subsBeha(:,2,1);
subs.s3_acc = subsBeha(:,3,1);

subs.s1_rt = subsBeha(:,1,2);
subs.s2_rt = subsBeha(:,2,2);
subs.s3_rt = subsBeha(:,3,2);

beha = subs;
save('beha.mat','beha')
%%

setStr = {'Load 1','Load 2','Load 4'};
indStr = {'Accuracy','RT','dprime'};
groupStr = {'Young','Old'};

clear X
X(:,1) = reshape(subsBeha(:,:,1),size(subsBeha,1)*size(subsBeha,2),1);
X(:,2) = repmat(subs.group,size(subsBeha,2),1);% between
X(:,3) = sort(repmat([1;2;3],subN,1));% within
X(:,4) = repmat([1:subN]',3,1);%subject
  [SSQs.acc, DFs.acc, MSQs.acc, Fs.acc, Ps.acc]=mixed_between_within_anova(X);
 
 clear X
X(:,1) = reshape(subsBeha(:,:,2),size(subsBeha,1)*size(subsBeha,2),1);
X(:,2) = repmat(subs.group,size(subsBeha,2),1);% between
X(:,3) = sort(repmat([1;2;3],subN,1));% within
X(:,4) = repmat([1:subN]',3,1);%subject
[SSQs.rt, DFs.rt, MSQs.rt, Fs.rt, Ps.rt]=mixed_between_within_anova(X);

%%
myFigBasic
myColors = flipud(crameri('bamako',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
mksize = 6;
figure('Name','Beha','Position',[200 200 400 400]);

t = 1;

subplot(2,2,t+2);hold all;axis square;
for g = 1:2
    tmpID = subs.group==g;
    dat= squeeze(subsBeha(tmpID,:,t));
    mn(g,:) = mean(dat);
    se(g,:) = std(dat)./sqrt(sum(tmpID));    
end

b = bar(mn','BarWidth',1);
b(1).FaceColor = myColors(1,:);
b(2).FaceColor = myColors(2,:);
for g = 1:2
    tmpID = subs.group==g;
    dat= squeeze(subsBeha(tmpID,:,t));
    scatter(b(g).XEndPoints,dat,mksize,'k','MarkerFaceColor',myColors(g,:))
end
errorbar(vertcat(b.XEndPoints)',mn',se','k.','LineWidth',1)

set(gca,'YLim',[0.6 1.09],'YTick',[0.6:0.2:1],'XTickLabel',setStr,'XLim',[0.25 3.75],'XTick',1:3)
legend(groupStr,'Location','southoutside')
ytickformat('%,.1f')
ylabel(indStr{t})

%---plot significance
[~,tmp_val] = ttest2(squeeze(subsBeha(subs.group==1,:,t)),squeeze(subsBeha(subs.group==2,:,t)));
sigY = [1 1 1]*1.03;
tmp_xcord = vertcat(b.XEndPoints)';

for i = 1:length(tmp_val)
    if tmp_val(i)<.05
        sigStr = '*';
        if tmp_val(i)<.01
            sigStr = '**';
            if tmp_val(i)<.001
                sigStr = '***';
            end
        end
        plot(tmp_xcord(i,:),[sigY(i) sigY(i)] ,'color','k','HandleVisibility','off')
        text(mean(tmp_xcord(i,:)),sigY(i),sigStr,'FontSize',18,'HorizontalAlignment','center')
    end
end


h4 = subplot(2,2,t);hold all;axis square
for g = 1:2
    tmpID = subs.group==g;
    dat = squeeze(subsBeha(tmpID,:,t));%
    boxplot(dat,'Positions',[2 6 8]+g*1.5,'Whisker',1,'Widths',0.8,'OutlierSize',3,'Symbol','ko')
end

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),myColors((j<=3)+1,:),'FaceAlpha',.7);
end

set(gca,'xtick',[2 6 10]+1.5*1.5,'XTickLabel',setStr)
ylabel(indStr{t})
set(gca,'YLim',[0.5 1],'XLim',[2 14.5])
ytickformat('%,.1f')
lines = findobj(h4, 'type', 'line');
set(lines, 'Color', 'k','linewidth',1.2);

t = 2;
subplot(2,2,t+2);hold all;axis square;
for g = 1:2
    tmpID = subs.group==g;
    dat= squeeze(subsBeha(tmpID,:,t));
    mn(g,:) = mean(dat);
    se(g,:) = std(dat)./sqrt(sum(tmpID));    
end
b = bar(mn','BarWidth',1);
b(1).FaceColor = myColors(1,:);
b(2).FaceColor = myColors(2,:);
for g = 1:2
    tmpID = subs.group==g;
    dat= squeeze(subsBeha(tmpID,:,t));
    scatter(b(g).XEndPoints,dat,mksize,'k','MarkerFaceColor',myColors(g,:))
end
errorbar(vertcat(b.XEndPoints)',mn',se','k.','LineWidth',1)
set(gca,'YLim',[0.4 1.6],'YTick',[0.4:0.4:3],'XTickLabel',setStr,'XLim',[0.5 3.5],'XTick',1:3)
legend(groupStr,'Location','southoutside')
ytickformat('%.1f')
ylabel(indStr{t})

%---plot significance
[~,tmp_val] = ttest2(squeeze(subsBeha(subs.group==1,:,t)),squeeze(subsBeha(subs.group==2,:,t)));
sigY = [1.4 1.4 1.4];
tmp_xcord = vertcat(b.XEndPoints)';

for i = 1:length(tmp_val)
    if tmp_val(i)<.05
        sigStr = '*';
        if tmp_val(i)<.01
            sigStr = '**';
            if tmp_val(i)<.001
                sigStr = '***';
            end
        end
        plot(tmp_xcord(i,:),[sigY(i) sigY(i)] ,'color','k','HandleVisibility','off')
        text(mean(tmp_xcord(i,:)),sigY(i),sigStr,'FontSize',18,'HorizontalAlignment','center')
    end
end

h4=subplot(2,2,t);hold all;axis square
dat = squeeze(subsBeha(:,:,t));%
for g = 1:2
    tmpID = subs.group==g;
    dat = squeeze(subsBeha(tmpID,:,t));%
    boxplot(dat,'Positions',[2 6 10]+g*1.5,'Whisker',1,'Widths',0.8,'OutlierSize',3,'Symbol','ko')
end

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),myColors((j<=3)+1,:),'FaceAlpha',.7);
end
set(gca,'xtick',[2 6 10]+1.5*1.5,'XTickLabel',setStr,'XLim',[2 14.5])
ylabel(indStr{t})
set(gca,'YLim',[0 2])
ytickformat('%,.1f')
lines = findobj(h4, 'type', 'line');
set(lines, 'Color', 'k','linewidth',1.2);


% t = 3;
% 
% subplot(2,3,t+3);hold all;axis square;box on
% for g = 1:2
%     tmpID = subs.group==g;
%     dat= squeeze(subsBeha(tmpID,:,t));
%     mn(g,:) = mean(dat);
%     se(g,:) = std(dat)./sqrt(sum(tmpID));    
% end
% b = bar(mn');
% b(1).FaceColor = myColors(1,:);
% b(2).FaceColor = myColors(2,:);
% for g = 1:2
%     tmpID = subs.group==g;
%     dat= squeeze(subsBeha(tmpID,:,t));
%     scatter(b(g).XEndPoints,dat,10,'k','MarkerFaceColor',myColors(g,:))
% end
% errorbar(vertcat(b.XEndPoints)',mn',se','k.','LineWidth',1)
% set(gca,'YLim',[1 5],'XTickLabel',setStr,'XLim',[0.5 3.5],'XTick',1:3)
% legend(groupStr,'Location','southoutside')
% ytickformat('%.1f')
% title(indStr{t})
% 
% %---plot significance
% [~,tmp_val] = ttest2(squeeze(subsBeha(subs.group==1,:,t)),squeeze(subsBeha(subs.group==2,:,t)));
% sigY = [4.7 4.7 4.7];
% tmp_xcord = vertcat(b.XEndPoints)';
% 
% for i = 1:length(tmp_val)
%     if tmp_val(i)<.05
%         sigStr = '*';
%         if tmp_val(i)<.01
%             sigStr = '**';
%             if tmp_val(i)<.001
%                 sigStr = '***';
%             end
%         end
%         plot(tmp_xcord(i,:),[sigY(i) sigY(i)] ,'color','k','HandleVisibility','off')
%         text(mean(tmp_xcord(i,:)),sigY(i),sigStr,'FontSize',18,'HorizontalAlignment','center')
%     end
% end
% 
% h4=subplot(2,3,t);hold all;axis square
% for g = 1:2
%     tmpID = subs.group==g;
%     dat = squeeze(subsBeha(tmpID,:,t));%
%     boxplot(dat,'Positions',[2 6 10]+g*1.5,'Whisker',1,'Widths',0.8,'OutlierSize',3,'Symbol','ko')
% end
% 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),myColors((j<=3)+1,:),'FaceAlpha',.7);
% end
% set(gca,'xtick',[2 6 10]+1.5*1.5,'XTickLabel',setStr)
% title(indStr{t})
% set(gca,'YLim',[0 6.5],'XLim',[2 14.5])
% ytickformat('%,.1f')
% lines = findobj(h4, 'type', 'line');
% set(lines, 'Color', 'k','linewidth',1.2);

saveas(gcf,fullfile(Dir.figs,'beha.png'))
saveas(gcf,fullfile(Dir.figs,'beha.pdf'))