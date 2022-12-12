clear;rng('shuffle');
load('subs.mat');

myFigBasic
addpath(genpath('D:\intWM-E\toolbox\fieldtrip-20211209'))
ft_defaults;

IsLap = 1;
IsdePhase = 1;
IsCorretTrials = 1;
txtCell = {'','','';'_lap','_dephase','_corrTrials'};
permFile = ['Perm_probePLV',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat'];
%% participants

tmp = find(subs.group==1 & ~ismember(subs.name,{'S05','S28','S27','SE08','SE16'}));

% S05: <20% EEG trials remained with a threshold of 75uV
% S27: quit after 4 blocks
% S28: bad data
mySubs{1} = tmp;% YOUNG

tmp = find(subs.group==2 & ~ismember(subs.name,{'S05','S28','S27','SE01','SE05','SE08','SE16'}));
mySubs{2} = tmp;% OLD
%%
clear tfAll
for gi = 1:2
    tmp = nan(length(mySubs{gi}),3);
    tmpB = nan(length(mySubs{gi}),3);
    parfor sub_i = 1:length(mySubs{gi})
        sn = mySubs{gi}(sub_i);
        subname = subs.name{sn};

        M = load(fullfile(Dir.results,[subname,'_thetaPLV',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']));

        PFCidx = contains(M.PLV.data{1}.label,'PFC');

        for cond_i = 1:3
            datAll = squeeze(mean(M.PLV.data{cond_i}.plvspctrm,4,'omitnan'));
            tmp(sub_i,cond_i) = mean(mean(datAll(PFCidx,~PFCidx)));% take the inter-regional PLV

            datAll = squeeze(mean(M.PLV.dataB{cond_i}.plvspctrm,4,'omitnan'));
            tmpB(sub_i,cond_i) = mean(mean(datAll(PFCidx,~PFCidx)));% take the inter-regional PLV

%             datAll = squeeze(mean(M.PLV.coh{cond_i}.cohspctrm,3));
%             tmp(sub_i,cond_i) = mean(mean(datAll(PFCidx,~PFCidx)));%
% 
%             datAll = squeeze(mean(M.PLV.cohB{cond_i}.cohspctrm,3));
%             tmpB(sub_i,cond_i) = mean(mean(datAll(PFCidx,~PFCidx)));% 

%             datAll = squeeze(mean(M.PLV.itpc{cond_i}.itpc,3));
%             tmp(sub_i,cond_i) = mean(mean(datAll(PFCidx,~PFCidx)));%
% 
%             datAll = squeeze(mean(M.PLV.itpc{cond_i}.itpcB,3));
%             tmpB(sub_i,cond_i) = mean(mean(datAll(PFCidx,~PFCidx)));% 

        end
    end
    tfAll.plv{gi} = tmp;
    tfAll.plvB{gi} = tmpB;
end

condStr = {'Fixation','Ignore','Attend','Attend+Resp'};
groupStr = {'Young','Older'};
%% PLV change

for gi = 1:2
    tfAll.plvDiff{gi} = tfAll.plv{gi};

    tfAll.plvDiff{gi} = tfAll.plv{gi}-tfAll.plvB{gi};
%     tfAll.plvDiff{gi} = tfAll.plv{gi}./tfAll.plvB{gi}-1;

end
% tfAll.plvDiff{1}([1],:) = [];

%% 2*3 mixed ANOVA



%%
myFigBasic
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
myColors = crameri('bam',4);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
figure;
axis square;hold all;
mn = cellfun(@mean,tfAll.plvDiff,'UniformOutput',false);
mn = vertcat(mn{:});

std = cellfun(@std,tfAll.plvDiff,'UniformOutput',false);
std = vertcat(std{:});
se = std./sqrt(cellfun(@length,mySubs'));

hb = bar(mn);
for i = 1:length(hb)
    hb(i).FaceColor = myColors(i,:);
end
xcoord = vertcat(hb.XEndPoints)';
errorbar(xcoord,mn,se,'k','LineStyle','none')
set(gca,'XTickLabel',groupStr,'xtick',[1 2])
legend(condStr)
ytickformat('%,.2f')

%%
% if rAnova.RotACC_Y.pValue(3)<.05
%     tmpA = linspace(0.75,0.85,6);
%     tmp_val = rotACC_tbl.Y.pValue([1 2 3 5 6 9]);
%     % fixation vs. ignore; fixation vs. attend; fixation-attend_res;
%     % ignore-attend; ignore-attend+res
%     tmp = vertcat(bar_h(:).XEndPoints);
%     tmp_xcord = {[1 2],[1 3],[1 4],[2:3],[2 4],3:4};
% 
%     for i = 1:length(tmp_val)
%         if tmp_val(i)<.05 
%             sigStr = '*';
%             if tmp_val(i)<.01
%                 sigStr = '**';
%                 if tmp_val(i)<.001
%                     sigStr = '***';
%                 end
%             end
%             sigY = max(max(dat))*tmpA(i);
%             plot(tmp_xcord{i}',[sigY sigY] ,'color','k','HandleVisibility','off')
%             text(mean(tmp_xcord{i}),sigY,sigStr,'FontSize',20,'HorizontalAlignment','center')
%         end
%     end
% end
%%

load Behav_Y.mat
wanted = ismember(str2num(demoInfo.name(:,2:end)),mySubs{1});
dat.beha{1} = squeeze(B.swapNT(wanted,:,1));% k of ignore and attend condition
dat.beha{1} = memResults.rotBegin(wanted,:);
dat.beha{1} = memResults.rotACC(wanted,:);
% dat.beha{1} = memResults.rotDif(wanted,1:3);

load Behav_E.mat
wanted = ismember(str2num(demoInfo.name(:,3:end)),mySubs{2}-28);
dat.beha{2} = squeeze(B.swapNT(wanted,:,1));% k of ignore and attend condition
dat.beha{2} = memResults.rotBegin(wanted,:);
dat.beha{2} = memResults.rotACC(wanted,:);
% dat.beha{2} = memResults.rotDif(wanted,:);

clear g
for gi = 1:2
    [r,pval] = corr(dat.beha{gi},tfAll.plvDiff{gi},'type','Spearman')
%     [r,pval] = corr(dat.beha{gi},dat.acc{gi},'type','Pearson')

    for condi = 1:3
        color_point=myColors(condi,:); %Set the color of the point,Three numbers respectively[R G B]Weight based on0~1Between
        X = tfAll.plvDiff{gi}(:,condi);
        Y = dat.beha{gi}(:,condi);
        g(gi,condi)=gramm('x',X,'y',Y); %Specify the values ??of the horizontal axis x and the vertical axis y, and create a gramm drawing object
        g(gi,condi).geom_point(); %Draw a scatter plot
        g(gi,condi).stat_glm(); %Draw lines and confidence intervals based on scatter plots
        g(gi,condi).set_names('x','PLV','y','Beha (k)'); %Set the title of the axis
        g(gi,condi).set_color_options('map',color_point); %Set the color of the point
        g(gi,condi).set_title(sprintf('%s-%s\nRho = %.3f, p = %.3f',groupStr{gi},condStr{condi},r(condi,condi),pval(condi,condi)));
        g(gi,condi).set_text_options('base_size' ,10,'title_scaling' ,1);%Set the font size, the base font size base_size is set16, the title font size of the axis is set to the base font size1.2Times
    end
end
figure('Position',[100 100 1200 500]);
g.draw();   
% saveas(gcf,fullfile('Figs',['',txtCell{IsCorretTrials+1,1},txtCell{IsOcci+1,2},txtCell{IsEOG+1,3},'.bmp']))


