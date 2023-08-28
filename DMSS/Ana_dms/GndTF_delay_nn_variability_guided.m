
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
varStr = {'LowVar','HighVar'};
% chanROI = {'CPz','CP1','CP2','Pz','Cz'};
chanROI = {'FCz','Fz','F1','F2','AFz'};

freq.betaFreq = [15 25];% Hz

timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s

%%

for Nsamps = 16%[16:4:40]
    clear subsAll
    tbl = table;

    for sub_i = 1:subN
        subname = subs.name{sub_i};
        disp(subname)
        matName =fullfile(Dir.results,[subname,'_delay_varGpow_NsampledAT',num2str(Nsamps),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
        if isfile(matName)
            load(matName);

            freqID.beta = dsearchn(tfDat{1}{1}.freq',freq.betaFreq');
            freqID.beta = freqID.beta(1):freqID.beta(2);

            timeID.all = dsearchn(tfDat{1}{1}.time',timeROI.all');
            timeID.all = timeID.all(1):timeID.all(2);

            timeID.Post = dsearchn(tfDat{1}{1}.time',timeROI.Post');
            timeID.Post = timeID.Post(1):timeID.Post(2);

            [log1,chanID] = ismember(chanROI,tfDat{1}{1}.label);

            for cond_i = 1:3
                for b = 1:length(tfDat{cond_i})
                    dat{sub_i}.betaAvgFrontalVar{cond_i}(b,1) = squeeze(mean(mean(mean(tfDat{cond_i}{b}.powspctrm(chanID(log1),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));
                end
                % Calculate the first and last quarters (25% and 75% percentiles)
                q1 = prctile(dat{sub_i}.betaAvgFrontalVar{cond_i}, 25);% low variability
                q4 = prctile(dat{sub_i}.betaAvgFrontalVar{cond_i}, 75);% high variability

                Quarter_Idx{1} = find(dat{sub_i}.betaAvgFrontalVar{cond_i} <= q1);
                Quarter_Idx{2} = find(dat{sub_i}.betaAvgFrontalVar{cond_i} >= q4);

                for q = 1:2
                    tmp = [];tmpB = [];
                    for i = 1:numel(Quarter_Idx{q})
                        tmp = [tmp;squeeze(mean(mean(mean(PowTrl{cond_i}{Quarter_Idx{q}(i)}.powspctrm(:,chanID(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"))];
                        tmpB = [tmpB;TrlBeha{cond_i}{Quarter_Idx{q}(i)}.button_resp_rt];
                    end
                    QuarteData.Power{sub_i,cond_i,q} = tmp;
                    QuarteData.Beha{sub_i,cond_i,q} = tmpB;

                    tmp_tbl = table(tmpPow,tmpBeha,'VariableNames',{'eeg','beha'});
                    tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);
                    tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
                    tmp_tbl.ss = ones(height(tmp_tbl),1)*cond_i;
                    tmp_tbl.quarter = ones(height(tmp_tbl),1)*q;

                    tbl = [tbl;tmp_tbl];

                    [subsAll.corrR(sub_i,cond_i,q), subsAll.corrp(sub_i,cond_i,q)] = corr(tmp,tmpB,'type','Spearman');
                end
            end
        end
    end

    save(fullfile(Dir.results,['DelayVarGuidedPow_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'tbl','dat','subsAll','-v7.3')
    %% glme
    load(fullfile(Dir.results,['DelayVarGuidedPow_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))
    %--------selected channels
    tbl.subj = categorical(tbl.subj);
    tbl.ss = categorical(tbl.ss);
    tbl.group = categorical(tbl.group);
    tbl.quarter = categorical(tbl.quarter);

    % first put everything as fixed effect, ss has no main effect or
    % interaction with other terms (cp chans)
    glme = fitglme(tbl,'beha ~ eeg*quarter*group*ss+(1|subj)');
    t = dataset2table(glme.anova);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[chanROI{:},'_anova1'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[chanROI{:},'_Coefficients1'])

    % control for ss, interaction group*eeg*quarter

    glme = fitglme(tbl,'beha ~ eeg*quarter*group+(1|subj)+(1|ss)');
    t = dataset2table(glme.anova);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[chanROI{:},'_anova2'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[chanROI{:},'_Coefficients2'])

    %----separate quater

    for q = 1:2
        tmp_tbl = tbl(double(tbl.quarter) == q,:) ;
        glme = fitglme(tmp_tbl,'beha ~ eeg*group+(1|subj)+(1|ss)');
        t = dataset2table(glme.anova);
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[varStr{q},'_anova'])
        t = dataset2table(glme.Coefficients);
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[varStr{q},'_Coefficients'])

        for gi = 1:2
            tmp_tbl = tbl(double(tbl.quarter) == 1 & double(tbl.group) == 2,:) ;
            glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)+(1|ss)');
            t = dataset2table(glme.anova);
            t.fomula{1} = char(glme.Formula);
            writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
                'FileType','spreadsheet','Sheet',[varStr{q},groupStr{gi},'_anova'])
            t = dataset2table(glme.Coefficients);
            t.fomula{1} = char(glme.Formula);
            writetable(t,fullfile(Dir.ana,'TableOutput',['DelayVariabilityGuidedPowerNNglme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
                'FileType','spreadsheet','Sheet',[varStr{q},groupStr{gi},'_Coefficients'])
        end
    end

    %% get the correlation coefficients of each subj and hist plot

    % histogram of correlation Rho
    myColors = flipud(crameri('bamako',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
    grougAligStr = {'Right','Left'};
    figure('position',[100 100 600 600])
    for q = 1:2
        for cond_i = 1:3
            subplot(2,3,cond_i+(q-1)*3);axis square;  hold all;

            for gi  = 1:2
                tmpsubs = find(subs.group == gi);
                [~,p,~,stats] = ttest(subsAll.corrR(tmpsubs,cond_i,q));

                h = histogram(subsAll.corrR(tmpsubs,cond_i,q));
                h.Normalization = 'probability';
                h.BinWidth = 0.14;
                h.FaceColor = myColors(gi,:);
                h.FaceAlpha = 0.9;
                plot([mean(subsAll.corrR(tmpsubs,cond_i,q)) mean(subsAll.corrR(tmpsubs,cond_i,q))],[0 1],'color',myColors(gi,:),'LineWidth',2,'HandleVisibility','off')
                text(2.2-1.1*gi,0.8,sprintf('Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.corrR(tmpsubs,cond_i,q)),stats.tstat,p,computeCohen_d(subsAll.corrR(tmpsubs,cond_i,q),zeros(length(tmpsubs),1),'paired')),...
                    'sc','color',myColors(gi,:),'HorizontalAlignment',grougAligStr{gi})
            end
            legend(groupStr,'Location','southoutside');
            ylabel('Proportion')
            xlabel('Spearman Rho value')
            title([chanROI{:}],condStr{cond_i})
            set(gca,'ylim',[0 0.5],'xlim',[-0.8 0.5])
        end
    end
    saveas(gca,fullfile(Dir.figs,['DelayVarianceGuidedPowerRTcorr_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
end