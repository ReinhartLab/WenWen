
clear;
load subsinfo.mat

RTsubs = nan(subN,2,6,6);
ACCsubs = RTsubs;
outPara.outRT = 6;% in s
outPara.outSD = 3; % SD

% R.slopePerBlock = nan(subN,2,6);


for sub_i = 1:subN
    if subs.exclude == 0
    subname = subs.name{sub_i};

    load(fullfile(Dir.beha.IBL,subname,[subname,'.mat']));
    subs.stim(sub_i) = subInfo.stim;
    subs.gender{sub_i} = subInfo.gender;
    subs.age(sub_i) = subInfo.age;
    subs.hand{sub_i} = subInfo.hand;
    ACC = nan(2,6,params.trlNperTaskperRun);% three stim phases * 6 blocks* repeat
    RT = ACC;
    clear eachTaskACC
    for p = 1:2 % part: 1 = pre, 2 = stim
        for bi = 1:6 % block
            tmp_table = subsMat(p).table(bi).all;
         
           
            for ri = 1:params.trlNperTaskperRun
                ACC(p,bi,ri) = mean(tmp_table.acc(tmp_table.repi == ri)==1);% averaged acc across three tasks
                if ACC(p,bi,ri)~=0

                    mn = mean(tmp_table.RT(tmp_table.repi == ri & tmp_table.acc == 1));
                    sd = std(tmp_table.RT(tmp_table.repi == ri & tmp_table.acc == 1));
                    outPara.up = mn+outPara.outSD*sd;
                    outPara.down = mn-outPara.outSD*sd;

                    tmpIdx = tmp_table.RT<outPara.outRT & tmp_table.RT < outPara.up & tmp_table.RT > outPara.down & tmp_table.repi == ri & tmp_table.acc == 1;
                    RT(p,bi,ri) = mean(tmp_table.RT(tmpIdx));
                end
            end

            for t = 1:3
                tmpIdx = tmp_table.curTask==t;
                ACCstat.eachTaskACC(sub_i,p,bi,t) = sum(tmp_table.acc(tmpIdx)==1)/params.trlNperTaskperRun;
            end

            x = 1:6;
            y = squeeze(RT(p,bi,:));
            tmp = polyfit(x(~isnan(y)),y(~isnan(y)),1);% remove nan
            R.slopePerBlock(sub_i,p,bi) = tmp(1);
        end

        x = 1:6;
        y = squeeze(mean(RT(p,:,:),2,'omitnan'));% array of trials,average across blocks and then fitting
        tmp_slope(sub_i,p,:) = polyfit(x,y,1);% slope and intercept
    end

    R.slopePre(sub_i) = tmp_slope(sub_i,1,1);
    R.slopeStim(sub_i) = tmp_slope(sub_i,2,1);

    RTsubs(sub_i,:,:,:) = RT;% subN*phase*block*trial
    ACCsubs(sub_i,:,:,:) = ACC;

    subs.faceACC(sub_i) = mean(ACCstat.eachTaskACC(sub_i,:,:,1),"all");
    subs.vehiACC(sub_i) = mean(ACCstat.eachTaskACC(sub_i,:,:,2),"all");
    subs.motionACC(sub_i) = mean(ACCstat.eachTaskACC(sub_i,:,:,3),"all");

    subs.faceACCpre(sub_i) = mean(ACCstat.eachTaskACC(sub_i,1,:,1),3);
    subs.vehiACCpre(sub_i) = mean(ACCstat.eachTaskACC(sub_i,1,:,2),3);
    subs.motionACCpre(sub_i) = mean(ACCstat.eachTaskACC(sub_i,1,:,3),3);

    R.PreRT(sub_i,:) = squeeze(mean(RT(1,:,:),2,'omitnan'))';% array of trials
    R.StimRT(sub_i,:) = squeeze(mean(RT(2,:,:),2,'omitnan'))';% array of trials
end
end
%%
R.slopePre= mean(R.slopePerBlock(:,1,:),3,'omitnan');
R.slopeStim= mean(R.slopePerBlock(:,2,:),3,'omitnan');

%% Demographic
fprintf(['\nfemale = %d,male = %d\n' ...
    'age:mean = %.1f, SD = %.1f' ...
    '\nhand: left = %d, right =%d\r'],sum(contains(subs.gender,'f')),sum(contains(subs.gender,'m')),mean(subs.age),std(subs.age),sum(contains(subs.hand,'l')),sum(contains(subs.hand,'r')));

%%
stimType = unique(subs.stim);

for gi = 1:2
    sub_idx = subs.stim == stimType(gi);
    RTstat.data{gi} = squeeze(mean(RTsubs(sub_idx,:,:,:),3,'omitnan'));%subj*part*trial
    RTstat.mn(gi,:,:) = mean(RTstat.data{gi});% cond*part*trial
    RTstat.se(gi,:,:) = std(RTstat.data{gi})./sqrt(sum(sub_idx));

    RTstat.bslData{gi} = bsxfun(@minus, RTstat.data{gi},RTstat.data{gi}(:,:,1));
    RTstat.bslmn(gi,:,:) = mean(RTstat.bslData{gi});% cond*part*trial
   
    ACCstat.data{gi} = mean(ACCsubs(sub_idx,:,:,:),3,'omitnan');
    ACCstat.mn(gi,:,:) = mean(ACCstat.data{gi});
    ACCstat.se(gi,:,:) = std(ACCstat.data{gi})./sqrt(sum(sub_idx));

    for bi = 1:6
        RTstat.blockdata{gi}{bi} = squeeze(RTsubs(sub_idx,:,bi,:));%subj*part*trial
        RTstat.blockmn(gi,bi,:,:) = mean(RTstat.blockdata{gi}{bi},'omitnan');% cond*block*part*trial
        RTstat.blockse(gi,bi,:,:) = std(RTstat.blockdata{gi}{bi},'omitnan')./sqrt(sum(sub_idx));
    end

    RTstat.slope{gi} = tmp_slope(sub_idx,:,:);% slope and intercept

end
subs.group = subs.stim;
subs.group(subs.group==2000) = 1;
subs.group(subs.group==2020) = 2;

% spssTable = subs(:,{'name','group','stim','slopePre','slopeStim','slopePost','PreRT','StimRT','PostRT'});
%%  plot pre-stim-post, average blocks
condStr = {'Sham','Theta','Beta'};

myColors = jet(6);
phaseStr = {'Pre','Stim','Post'};
blockStr = {'Block1','Block2','Block3','Block4','Block5','Block6'};

figure('Position',[0 0 1300 1000]);

for gi = 1:3
    for p = 1:3
        subplot(3,3,p);hold all;

        plot(squeeze(RTstat.mn(gi,p,:)),'-','Color',myColors(gi,:),'LineWidth',3)
        errorbar(squeeze(RTstat.mn(gi,p,:)),squeeze(RTstat.se(gi,p,:)),'Color',myColors(gi,:),'LineWidth',3,'HandleVisibility','off')
        %         for j = 1:sum(sub_idx)
        %             %             plot(squeeze(RTstat.data{cond_i}(j,si,:)),'--','Color',myColors(cond_i,:),'HandleVisibility','off')
        %             plot(squeeze(RTstat.data{cond_i}(j,si,:)),'--','HandleVisibility','off')
        %         end

        plot(squeeze(ACCstat.mn(gi,p,:)),'o--','Color',myColors(gi,:),'LineWidth',3,'HandleVisibility','off')
        xlabel('Trial')
        ylabel('Accuracy&RT(s)');
        title(phaseStr{p})
        set(gca,'YLim',[0.7 2],'YTick',[0.5:0.5:3])


        subplot(3,3,p+3);hold all;
        errorbar(squeeze(RTstat.bslmn(gi,p,:)),squeeze(RTstat.se(gi,p,:)),'Color',myColors(gi,:),'LineWidth',3)
        %         for j = 1:sum(sub_idx)
        %             %             plot(squeeze(RTstat.data{cond_i}(j,si,:)),'--','Color',myColors(cond_i,:),'HandleVisibility','off')
        %             plot(squeeze(RTstat.data{cond_i}(j,si,:)),'--','HandleVisibility','off')
        %         end
        legend(condStr,'Location','northoutside')

        xlabel('Trial')
        ylabel('RT changes(s)');
        title(phaseStr{p},'Normalized')
        set(gca,'YLim',[-0.5 0.2],'YTick',[-0.6:0.2:3])

    end
end
for gi = 1:3
    for p = 1:3
        subplot(3,3,p+6);hold all;

        errorbar(squeeze(RTstat.rMean(gi,p,:)),squeeze(RTstat.rSE(gi,p,:)),'Color',myColors(gi,:),'LineWidth',3,'HandleVisibility','off')
        xlabel('Trial')
        ylabel('Running RT(s)');
        title(phaseStr{p},'Running RT')
        set(gca,'YLim',[1 2],'YTick',[0.5:0.5:3])
    end
end
saveas(gcf,fullfile('FiguresOutput','Beha_RTacc.png'))
%% plot: pre-stim-post each block
figure('Position',[50 0 1500 800]);

for gi = 1:3
    for p = 1:3
        subplot(3,3,p+(gi-1)*3);
        hold all;

        for bi = 1:6
            errorbar(squeeze(RTstat.blockmn(gi,bi,p,:)),squeeze(RTstat.blockse(gi,bi,p,:)),'.-','Color',myColors(bi,:),'LineWidth',3)
        end
        xlabel('Trial')
        ylabel('RT(s)');
        title([phaseStr{p},condStr{gi}],sprintf('N=%d',sum(subs.stim == stimType(gi))))
        set(gca,'YLim',[0.7 2.2],'YTick',[0.5:0.5:2.5],'xlim',[0 7])
        legend(blockStr,'Location','eastoutside')
    end
end
saveas(gcf,fullfile('FiguresOutput','Beha_RTblockwise.png'))

% %% mixed anova on slope 3(pre-stim-post)*3(beta-theta-sham)
%
% clear X
% X(:,1) = [R.slopePre;R.slopeStim;R.slopePost];
% X(:,2) = [subs.group;subs.group;subs.group];
% X(:,3) = [ones(subN,1);ones(subN,1)*2;ones(subN,1)*3];
% X(:,4) = [1:subN 1:subN 1:subN]';
%
% [SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(X);

%% fitrm mixed anova


t = subs(:,{'group','slopePre','slopeStim','slopePost'});
t.group = categorical(t.group);
Meas = table([1 2 3]','VariableNames',{'Phase'});
rm = fitrm(t,'slopePre-slopePost~group','WithinDesign',Meas);
ranova(rm)% main within & interaction
anova(rm) % main between
% manova(rm,'By','group')
% manova(rm)

% corrctionStr = 'bonferroni';
corrctionStr = 'hsd';%tukey-kramer
% corrctionStr = 'lsd';%no correction


tbl = multcompare(rm,'Phase','ComparisonType',corrctionStr);% if main effect
tbl = multcompare(rm,'group','ComparisonType',corrctionStr);% if main effect
T_phase = multcompare(rm,'Phase','By','group','ComparisonType',corrctionStr); % interaction
T_phase.group=double(T_phase.group);

T_phase.cohenD(3) = computeCohen_d(R.slopeStim(subs.group==T_phase.group(3)),R.slopePre(subs.group==T_phase.group(3)),'paired');
T_phase.cohenD(9) = computeCohen_d(R.slopeStim(subs.group==T_phase.group(9)),R.slopePre(subs.group==T_phase.group(9)),'paired');
T_phase.cohenD(15) = computeCohen_d(R.slopeStim(subs.group==T_phase.group(15)),R.slopePre(subs.group==T_phase.group(15)),'paired');

T_phase.cohenD(4) = computeCohen_d(R.slopeStim(subs.group==T_phase.group(4)),R.slopePost(subs.group==T_phase.group(4)));
T_phase.cohenD(10) = computeCohen_d(R.slopeStim(subs.group==T_phase.group(10)),R.slopePost(subs.group==T_phase.group(10)));
T_phase.cohenD(16) = computeCohen_d(R.slopeStim(subs.group==T_phase.group(16)),R.slopePost(subs.group==T_phase.group(16)));

T_group = multcompare(rm,'group','By','Phase','ComparisonType',corrctionStr);
T_group.group_1 = double(T_group.group_1);
T_group.group_2 = double(T_group.group_2);

T_group.cohenD(8) = computeCohen_d(R.slopeStim(subs.group==T_group.group_1(8)),R.slopeStim(subs.group==T_group.group_2(8)),'independent');
T_group.cohenD(10) = computeCohen_d(R.slopeStim(subs.group==T_group.group_1(10)),R.slopeStim(subs.group==T_group.group_2(10)),'independent');

%% plot: learning slope

plotStat = 1;% whether to plot stats in figure;

figure('position',[200 200 600 600]);
hold all;axis square
for gi = 1:3
    plot(squeeze(RTstat.slope{gi}(:,:,1))','--','Color',myColors(gi,:),'HandleVisibility','off')% individuals
end
for gi = 1:3
    mn = squeeze(mean(RTstat.slope{gi}(:,:,1)));
    se =  squeeze(std(RTstat.slope{gi}(:,:,1)))./sqrt(size(RTstat.slope{gi},1));
    errorbar(mn,se,'.-','Color',myColors(gi,:),'LineWidth',5)

    text(0.05,0.8+gi*0.05,sprintf('%dHz N=%d',stimType(gi),size((RTstat.slope{gi}),1)),'sc')
end
ylabel('Learning slope');
set(gca,'XLim',[0 4],'xTick',1:3,'XTickLabel',phaseStr)
legend(condStr,'Location','eastoutside')
title(corrctionStr)
if plotStat
    % sham vs beta  \n  theta vs beta
    text(0.4,0.9,sprintf('Stim\nBeta vs Sham:p=%.3f,d = %.3f\nBeta vs theta:p=%.3f,d = %.3f',T_group.pValue(8),T_group.cohenD(8),T_group.pValue(10),T_group.cohenD(10)),'sc')% stim

    % sham: Pre vs Stim ; Post vs Stim
    text(0.05,0.15,sprintf('Sham: Pre-Stim p=%.3f,d = %.3f; Post-Stim p=%.3f,d = %.3f',T_phase.pValue(3),T_phase.cohenD(3),T_phase.pValue(4),T_phase.cohenD(4)),'sc')
    text(0.05,0.1,sprintf('Theta: Pre-Stim p=%.3f,d = %.3f; Post-Stim p=%.3f,d = %.3f',T_phase.pValue(9),T_phase.cohenD(9),T_phase.pValue(10),T_phase.cohenD(10)),'sc')
    text(0.05,0.05,sprintf('Beta: Pre-Stim p=%.3f,d = %.3f; Post-Stim p=%.3f,d = %.3f',T_phase.pValue(15),T_phase.cohenD(15),T_phase.pValue(16),T_phase.cohenD(16)),'sc')
end
% saveas(gcf,fullfile('FiguresOutput',['Beha_LearningSlope',corrctionStr '.pdf']))
saveas(gcf,fullfile('FiguresOutput',['Beha_LearningSlope',corrctionStr,'.png']))

%% baseline slope and baseline RT
mdl = fitlm(mean(R.PreRT,2),R.slopePre);
anova(mdl,'summary')
figure('Position',[30 30 450 250]);hold all
plot(mdl,'linewidth',2);
text(0.6,0.1,sprintf('p=%.3f',mdl.anova.pValue(1)),'sc')% stim
legend('location','bestoutside')
xlabel('Pre RT')
ylabel('Pre slope')
set(gca,'XLim',[0.5 3],'YLim',[-0.2 0.05])
title('RT & Slope',['N=',num2str(subN)])
saveas(gcf,fullfile('FiguresOutput',['Beha_PreSlopeCorrPreRT.png']))

%% stimulation slope and baseline slope

figure('Position',[30 30 1200 350]);
for cond_i = 1:3

    mdl = fitlm(R.slopePre(subs.group==cond_i),R.slopeStim(subs.group==cond_i));
    anova(mdl,'summary')

    ax1 = subplot(1,3,cond_i);
    plot(ax1,mdl,'linewidth',2);
    hold all
    text(0.6,0.1,sprintf('p=%.3f',mdl.anova.pValue(1)),'sc')% stim
    legend('location','southoutside')
    xlabel('Pre slope')
    ylabel('Stim slope')
    set(ax1,'XLim',[-0.18 0.06],'YLim',[-0.18 0.15],'XTick',-0.18:0.06:0.06,'YTick',-0.15:0.15:0.15)
    title(condStr{cond_i})
end
saveas(gcf,fullfile('FiguresOutput',['Beha_stimSlopeCorrPreSlope.png']))

%% glme
tbl = array2table(X,"VariableNames",{'slope','group','phase','subject'});
tbl.group = tbl.group; % reverse sham and beta
tbl.group = categorical(tbl.group);
tbl.phase = categorical(tbl.phase);
tbl.subject=categorical(tbl.subject);
glme = fitglme(tbl,'slope~group*phase+(1|subject)','DummyVarCoding','effects')

%%
tbl = table;
tbl.slope = reshape(R.slopePerBlock,subN*3*6,1);%subN-pre-block1 subN-stim-block1
tbl.group = repmat(subs.group,3*6,1);
tbl.group = categorical(tbl.group);
tbl.phase = repmat([ones(subN,1)*1;ones(subN,1)*2;ones(subN,1)*3],6,1);
tbl.phase = categorical(tbl.phase);
tbl.block = sort(repmat([1:6]',subN*3,1));
tbl.subject = repmat([1:subN]',3*6,1);
tbl.subject=categorical(tbl.subject);

% glme = fitglme(tbl,'slope~group*phase+(1+block|subject)')
% glme = fitglme(tbl,'slope~group*phase+(1|subject)','DummyVarCoding','full')
glme = fitglme(tbl,'slope~phase*group+(1|subject)','DummyVarCoding','effects')
