clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
addpath(genpath('D:\intWM-E\toolbox\gramm-master'))
addpath(genpath('D:\intWM-E\toolbox\crameri_v1.08'))

IsCorretTrials = 1;
IsBL2preDelay = 1;% baselined to pre response
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Both'};
condStr = {'ss1','ss2','ss4'};

% frontalROI = {'CPz','CP1','CP2','Pz','Cz'};
frontalROI = {'FCz','Fz','F1','F2','AFz'};

% occipROI = {'POz','Oz'};
freq.betaFreq = [15 25];% Hz
% freq.alphaFreq = [8 13];% Hz

timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s
timeROI.Pre = [-0.4 0];% in s
permN = 10000;
%%
for Nsamps = 40%[16:4:40]
    clear subsAll
    for sub_i = 1:subN
        subname = subs.name{sub_i};
        matName =fullfile(Dir.results,[subname,'_resp_ft_power_NsampledAT',num2str(Nsamps),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

        if isfile(matName)
            load(matName);

            freqID.beta = dsearchn(tfDat{1}{1}.freq',freq.betaFreq');
            freqID.beta = freqID.beta(1):freqID.beta(2);

            timeID.all = dsearchn(tfDat{1}{1}.time',timeROI.all');
            timeID.all = timeID.all(1):timeID.all(2);

            timeID.Post = dsearchn(tfDat{1}{1}.time',timeROI.Post');
            timeID.Post = timeID.Post(1):timeID.Post(2);

            timeID.Pre = dsearchn(tfDat{1}{1}.time',timeROI.Pre');
            timeID.Pre = timeID.Pre(1):timeID.Pre(2);

            clear chanID
            [log1,chanID.frontal] = ismember(frontalROI,tfDat{1}{1}.label);

            for cond_i = 1:3
                for b = 1:length(tfDat{cond_i})
                    dat{sub_i}.betaAvgFrontal{cond_i}(b,1) = squeeze(mean(mean(mean(tfDat{cond_i}{b}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));
                    dat{sub_i}.betaChans{cond_i}(b,:) = squeeze(mean(mean(tfDat{cond_i}{b}.powspctrm(:,freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"));

                    dat{sub_i}.behaACC{cond_i}(b,1) = mean(PostTrlBeha{cond_i}{b}.button_resp_corr,'omitnan');
                    dat{sub_i}.behaRT{cond_i}(b,1) = mean(PostTrlBeha{cond_i}{b}.button_resp_rt,'omitnan');
                end

                subsAll.betaFrontalVar(sub_i,cond_i) = mean(dat{sub_i}.betaAvgFrontal{cond_i},'omitnan');
                subsAll.behaACC(sub_i,cond_i) = mean(dat{sub_i}.behaACC{cond_i},'omitnan');
            end
        end
    end
    save(fullfile(Dir.results,['RespBetaPower_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll','dat','-v7.3')
    load(fullfile(Dir.results,['RespBetaPower_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))

    %% regress model collapse groups: beha(ACC) = betaVar*set size+(1|subj)+ (1|group)

    myfun = @(x) [x.behaACC{1};x.behaACC{2};x.behaACC{3}];
    accArray =  cellfun(myfun,dat,'UniformOutput',false);
    accArray =  vertcat(accArray{:});% sub1(ss1-ss2-ss4)-sub2(ss1-ss2-ss4)...

    for c = 1:size(dat{1}.betaChans{1},2)
        loadArray = [];
        subArray = [];
        betaArray = [];
        groupArray = [];

        for sub_i = 1:subN
            load1 = dat{sub_i}.betaChans{1}(:,c);
            load2 = dat{sub_i}.betaChans{2}(:,c);
            load3 = dat{sub_i}.betaChans{3}(:,c);
            tmp = [load1;load2;load3];

            betaArray = [betaArray;tmp];
            subArray = [subArray;ones(length(tmp),1)*sub_i];
            loadArray = [loadArray;ones(length(load1),1)*1;ones(length(load2),1)*2;ones(length(load3),1)*3];
            groupArray = [groupArray;ones(length(tmp),1)*subs.group(sub_i)];
        end
        % glme

        tbl = table(accArray,betaArray, subArray, loadArray,groupArray,'VariableNames',{'beha','eeg','subj','ss','group'});
        tbl.subj = categorical(tbl.subj);
        tbl.group = categorical(tbl.group);
        glme = fitglme(tbl,'beha ~ 1+eeg*ss+(1|subj)+(1|group)');

        t = dataset2table(glme.anova);
        myReg.groupBoth.glme.pBeta(c) = t.pValue(2);
        myReg.groupBoth.glme.pInter(c) = t.pValue(4);

        t = dataset2table(glme.Coefficients);
        myReg.groupBoth.glme.FBeta(c) = t.Estimate(2);
        myReg.groupBoth.glme.FInter(c) = t.Estimate(4);
    end

    %% separate for young and old

    for gi = 1:2
        tmpsubs = find(subs.group == gi);
        myfun = @(x) [x.behaACC{1};x.behaACC{2};x.behaACC{3}];
        accArray =  cellfun(myfun,dat(tmpsubs),'UniformOutput',false);
        accArray =  vertcat(accArray{:});% sub1(ss1-ss2-ss4)-sub2(ss1-ss2-ss4)...

        for c = 1:size(dat{1}.betaChans{1},2)
            loadArray = [];
            subArray = [];
            betaArray = [];

            for sub_i = 1:length(tmpsubs)
                load1 = dat{tmpsubs(sub_i)}.betaChans{1}(:,c);
                load2 = dat{tmpsubs(sub_i)}.betaChans{2}(:,c);
                load3 = dat{tmpsubs(sub_i)}.betaChans{3}(:,c);
                tmp = [load1;load2;load3];

                betaArray = [betaArray;tmp];
                subArray = [subArray;ones(length(tmp),1)*sub_i];
                loadArray = [loadArray;ones(length(load1),1)*1;ones(length(load2),1)*2;ones(length(load3),1)*3];
            end

            tbl = table(accArray,betaArray, subArray, loadArray,'VariableNames',{'beha','eeg','subj','ss'});
            tbl.subj = categorical(tbl.subj);
            glme = fitglme(tbl,'beha ~ 1+eeg*ss+(1|subj)');
            t = dataset2table(glme.anova);
            myReg.group{gi}.glme.pBeta(c) = t.pValue(2);

            t = dataset2table(glme.Coefficients);
            myReg.group{gi}.glme.FBeta(c) = t.Estimate(2);
        end
    end

    %%   plot B of beta
    myFigBasic;

    load chanLocs_dms.mat
    [is, loc] = ismember(tfDat{1}{1}.label,{chanloc.labels});
    chanloc  = chanloc(loc(is));
    threshP = 0.001;
    mkSize = 5;
%     myColors = jet;
        myColors = hot;
    %     myColors = myColors(256/2:256,:);
    % myColors = crameri('bam',256);

    figure('Position',[100 20 1000 500]);

    for gi = 1:2
        subplot(1,4,gi)
        topoplot(myReg.group{gi}.glme.FBeta,chanloc,'emarker2',{find(myReg.group{gi}.glme.pBeta<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
        % caxis([-5, 5])
        h = colorbar;
        h.Label.String = sprintf('GLME Coef(p<%.3f)',threshP);
        h.Label.Rotation = -90;
        h.Label.Position = h.Label.Position + [1 0 0];
        title(groupStr{gi},['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
    end
    subplot(1,4,3)
    topoplot(myReg.groupBoth.glme.FBeta,chanloc,'emarker2',{find(myReg.groupBoth.glme.pBeta<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    % caxis([-5, 5])
    h = colorbar;
    h.Label.String = sprintf('GLME Coef(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title('Collapse groups: eeg',['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);

    subplot(1,4,4)
    topoplot(myReg.groupBoth.glme.FInter,chanloc,'emarker2',{find(myReg.groupBoth.glme.pInter<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
    % caxis([-5, 5])
    h = colorbar;
    h.Label.String = sprintf('GLME Coef(p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1 0 0];
    title('Collapse groups: eeg*ss',['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
    saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoGLME_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))

    %% use ft to plot, so that pdf can be editable
    ft_temp=load(fullfile(Dir.results,['ERP',txtCell{1,1},txtCell{1,2},'.mat']));

    cfg= [];
    cfg.latency = 0;
    cfg.channel = tfDat{1}{1}.label;
    ft_temp = ft_timelockgrandaverage(cfg,ft_temp.subsAll{1:2,1});

    ft_temp.avg = myReg.groupBoth.glme.FBeta;
    %%
    figure('position',[100 100 350 350])
    cfg = [];
     cfg.parameter = 'avg';
     cfg.colormap = hot;
%        cfg.colormap = jet;
  cfg.highlight          =  'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize   = 20;
%     cfg.highlightchannel = tfDat{1}{1}.label(myReg.groupBoth.glme.pBeta<threshP & myReg.groupBoth.glme.FBeta>0);
        cfg.highlightchannel = tfDat{1}{1}.label(myReg.groupBoth.glme.pBeta<threshP);
% cfg.style = 'straight';
cfg.style = 'both';
    cfg.marker = 'off';
    cfg.layout = 'easycapM11.mat';
    cfg.comment = 'no';
    cfg.figure = 'gca';
    ft_topoplotER(cfg,ft_temp);
    %
    %     topoplot(myReg.groupBoth.glme.FBeta,chanloc,'emarker2',{find(myReg.groupBoth.glme.pBeta<threshP),'*','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',flipud(hot),'style','map');
    caxis([-1*10^-6, 1*10^-6])
    h = colorbar;
    h.Label.String = sprintf('GLME Coefficient (p<%.3f)',threshP);
    h.Label.Rotation = -90;
    h.Label.Position = h.Label.Position + [1.2 0 0];
    h.Ticks = [-1*10^-6 0 1*10^-6];
    h.TickLength = 0;
    h.FontSize = 13;
    title('N+1 correlation (Pow)',['ClusterSize=',num2str(Nsamps)],'FontSize',13);

    saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoGLMEn+1Coef_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.pdf']))

    %% separate groups and set size
    figure('Position',[100 20 1000 800]);
    for gi = 1:2
        tmpsubs = find(subs.group == gi);
        for s = 1:3
            for c = 1:size(dat{1}.betaChans{1},2)
                subArray = [];
                betaArray = [];
                accArray = [];

                for sub_i = 1:length(tmpsubs)
                    tmp = dat{tmpsubs(sub_i)}.betaChans{s}(:,c);
                    accArray = [accArray;dat{tmpsubs(sub_i)}.behaACC{s}];
                    betaArray = [betaArray;tmp];
                    subArray = [subArray;ones(length(tmp),1)*tmpsubs(sub_i)];
                end

                tbl = table(accArray,betaArray, subArray,'VariableNames',{'beha','eeg','subj'});
                tbl.subj = categorical(tbl.subj);
                glme = fitglme(tbl,'beha ~ eeg+(1|subj)');% random intercept with a fixed slope
                t = dataset2table(glme.anova);
                myReg.ss{gi,s}.glme.pBeta(c) = t.pValue(2);

                t = dataset2table(glme.Coefficients);
                myReg.ss{gi,s}.glme.FBeta(c) = t.Estimate(2);
            end
            subplot(2,3,3*(gi-1)+s);
            topoplot(myReg.ss{gi,s}.glme.FBeta,chanloc,'emarker2',{find(myReg.ss{gi,s}.glme.pBeta<threshP),'o','k',mkSize},'plotrad',0.53,'electrodes','off','colormap',myColors );
            % caxis([-5, 5])
            h = colorbar;
            h.Label.String = sprintf('GLME Coeff(p<%.3f)',threshP);
            h.Label.Rotation = -90;
            h.Label.Position = h.Label.Position + [1 0 0];
            title([groupStr{gi} condStr{s}],['N+1 corr ACC TrlCluster=',num2str(Nsamps)]);
        end
    end
    saveas(gcf,fullfile(Dir.figs,['RespBetaPowerTopoGLMEeachSS_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4} '.png']))
    %% examine channels that showed beta rebound post response using glme
    % Y(n+1) = set size*eeg(n)+(1|subj)+(1|group)
    % collapse groups
    myfun = @(x) [x.behaACC{1};x.behaACC{2};x.behaACC{3}];
    accArray =  cellfun(myfun,dat,'UniformOutput',false);
    accArray =  vertcat(accArray{:});% sub1(ss1-ss2-ss4)-sub2(ss1-ss2-ss4)...

    loadArray = [];
    subArray = [];
    betaArray = [];
    groupArray = [];
    for sub_i = 1:subN
        load1 = dat{sub_i}.betaAvgFrontal{1};
        load2 = dat{sub_i}.betaAvgFrontal{2};
        load3 = dat{sub_i}.betaAvgFrontal{3};
        tmp = [load1;load2;load3];

        betaArray = [betaArray;tmp];
        subArray = [subArray;ones(length(tmp),1)*sub_i];
        loadArray = [loadArray;ones(length(load1),1)*1;ones(length(load2),1)*2;ones(length(load3),1)*3];
        groupArray = [groupArray;ones(length(tmp),1)*subs.group(sub_i)];
    end

    tbl = table(accArray,betaArray, subArray, loadArray,groupArray,'VariableNames',{'beha','eeg','subj','ss','group'});
    tbl.subj = categorical(tbl.subj);
    tbl.group = categorical(tbl.group);
    myReg.groupBothFrontal.glme = fitglme(tbl,'beha ~ 1+eeg*ss+(1|subj)+(1|group)');

    %  t = dataset2table(myReg.groupBothFrontal.glme.anova)
    %% selected channels each set size
    threshP = 0.01;

    for s = 1:3
        subArray = [];
        betaArray = [];
        accArray = [];
        groupArray = [];

        for sub_i = 1:subN
            tmp = dat{sub_i}.betaAvgFrontal{s};
            accArray = [accArray;dat{sub_i}.behaACC{s}];
            betaArray = [betaArray;tmp];
            subArray = [subArray;ones(length(tmp),1)*sub_i];
            groupArray = [groupArray;ones(length(tmp),1)*subs.group(sub_i)];
        end

        tbl = table(accArray,betaArray, subArray,groupArray,'VariableNames',{'beha','eeg','subj','group'});
        tbl.subj = categorical(tbl.subj);
        tbl.group = categorical(tbl.group);
        myReg.ssFrontal{s}.glme = fitglme(tbl,'beha ~ eeg+(1|subj)+(1|group)');% random intercept with a fixed slope
        t = dataset2table(myReg.ssFrontal{s}.glme.anova);
        %     myReg.ss{s}.glme.pBeta= t.pValue(2);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),'FileType','spreadsheet','Sheet',[condStr{s},'_anova'])
        t = dataset2table(myReg.ssFrontal{s}.glme.anova);
        %
        t = dataset2table(myReg.ssFrontal{s}.glme.Coefficients);
        %     myReg.ss{s}.glme.FBeta = t.Estimate(2);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespBetaPowerN+1glme_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),'FileType','spreadsheet','Sheet',[condStr{s},'_Coefficients'])
    end
end