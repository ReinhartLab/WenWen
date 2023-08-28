
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

for Nsamps = 40
    clear subsAll
    tbl = table;

    %    tmpSubj = [];
    %     tmpPow = [];
    %     tmpSSn = []; % ss of n
    %     tmpBeha = [];
    %     tmpGroup = [];
    %     tmpQ = [];
    %     tmpSSn1 = []; % ss of n+1
    for sub_i = 1:subN
        subname = subs.name{sub_i};
        disp(subname)
        matName =fullfile(Dir.results,[subname,'_resp_varGpow_seq_NsampledAT',num2str(Nsamps),txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

        if isfile(matName)
            load(matName);

            freqID.beta = dsearchn(tfDat{1}.freq',freq.betaFreq');
            freqID.beta = freqID.beta(1):freqID.beta(2);

            timeID.all = dsearchn(tfDat{1}.time',timeROI.all');
            timeID.all = timeID.all(1):timeID.all(2);

            timeID.Post = dsearchn(tfDat{1}.time',timeROI.Post');
            timeID.Post = timeID.Post(1):timeID.Post(2);

            [log1,chanID] = ismember(chanROI,tfDat{1}.label);

            for b = 1:length(tfDat)
                dat{sub_i}.betaAvgFrontalVar(b,1) = squeeze(mean(mean(mean(tfDat{b}.powspctrm(chanID(log1),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));
            end
            % Calculate the first and last quarters (25% and 75% percentiles)
            q1 = prctile(dat{sub_i}.betaAvgFrontalVar, 25);% low variability
            q4 = prctile(dat{sub_i}.betaAvgFrontalVar, 75);% high variability

            Quarter_Idx{1} = find(dat{sub_i}.betaAvgFrontalVar <= q1);
            Quarter_Idx{2} = find(dat{sub_i}.betaAvgFrontalVar >= q4);

            for q = 1:2
                tmp = [];tmpB = [];tmpS1 = [];tmpS0 = [];
                for i = 1:numel(Quarter_Idx{q})
                    tmp = [tmp;squeeze(mean(mean(mean(PowTrl{Quarter_Idx{q}(i)}.powspctrm(:,chanID(log1),freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"),4,"omitnan"))];
                    tmpB = [tmpB;PostTrlBeha{Quarter_Idx{q}(i)}.button_resp_rt];
                    tmpS1 = [tmpS1;PostTrlBeha{Quarter_Idx{q}(i)}.ss_num];
                    tmpS0 = [tmpS0;TrlBeha{Quarter_Idx{q}(i)}.ss_num];
                end
                QuarteData.Power{sub_i,q} = tmp;
                QuarteData.Beha{sub_i,q} = tmpB;

                tmp_tbl = table(tmpPow,tmpBeha,tmpS1,tmpS0,'VariableNames',{'eeg','beha','ssn1','ssn'});
                tmp_tbl.group = ones(height(tmp_tbl),1)*subs.group(sub_i);
                tmp_tbl.subj = ones(height(tmp_tbl),1)*sub_i;
                tmp_tbl.quarter = ones(height(tmp_tbl),1)*q;

                tbl = [tbl;tmp_Tbl];

                [subsAll.corrR(sub_i,q), subsAll.corrp(sub_i,q)] = corr(tmp,tmpB,'type','Spearman');
            end
        end
    end
    save(fullfile(Dir.results,['RespVarGuidedPow_seqRT_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'tbl','dat','subsAll','-v7.3')
    %% glme
    load(fullfile(Dir.results,['RespVarGuidedPow_seqRT_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']))
    %--------selected channels
    tbl.subj = categorical(tbl.subj);
    tbl.ssn = categorical(tbl.ssn);
    tbl.ssn1 = categorical(tbl.ssn1);
    tbl.group = categorical(tbl.group);
    tbl.quarter = categorical(tbl.quarter);

    % first put everything as fixed effect

    glme = fitglme(tbl,'beha ~ eeg*quarter*group*ssn*ssn1+(1|subj)');
    t = dataset2table(glme.anova)
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[chanROI{:},'_anova1'])
    t = dataset2table(glme.Coefficients);
    t.fomula{1} = char(glme.Formula);
    writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        'FileType','spreadsheet','Sheet',[chanROI{:},'_Coefficients1'])

    %% ----separate quarter CP

    for q = 1:2
        tmp_tbl = tbl(double(tbl.quarter) == q,:) ;
        %         glme = fitglme(tmp_tbl,'beha ~ eeg*group*ssn*ssn1+(1|subj)');
        %                 glme = fitglme(tmp_tbl,'beha ~ eeg*group+(1|subj)');
        glme = fitglme(tmp_tbl,'beha ~ eeg*group+(1|subj)+(1|ssn1)');

        %                 glme = fitglme(tmp_tbl,'beha ~ eeg*group+(1|subj)+(1|ssn)+(1|ssn1)');
        t = dataset2table(glme.anova);
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[varStr{q},'_anova'])
        t = dataset2table(glme.Coefficients);
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[varStr{q},'_Coefficients'])

        %        for s = 1:3
        %             tmp_tbl = tbl(double(tbl.quarter) == q & double(tbl.ssn1) == s,:) ;
        %             glme = fitglme(tmp_tbl,'beha ~ eeg*ssn*group+(1|subj)');
        %
        %             t = dataset2table(glme.anova);
        %             t.fomula{1} = char(glme.Formula);
        %             writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        %                 'FileType','spreadsheet','Sheet',[varStr{q},condStr{s},'_anova'])
        %             t = dataset2table(glme.Coefficients);
        %             t.fomula{1} = char(glme.Formula);
        %             writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
        %                 'FileType','spreadsheet','Sheet',[varStr{q},condStr{s},'_Coefficients'])
        %        end

        for gi = 1:2
            tmp_tbl = tbl(double(tbl.quarter) == q & double(tbl.group) == gi,:) ;
            glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)');
            %               glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)+(1|ssn1)');
            %                 glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)+(1|ssn)+(1|ssn1)');
            t = dataset2table(glme.anova);
            t.fomula{1} = char(glme.Formula);
            writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
                'FileType','spreadsheet','Sheet',[varStr{q},groupStr{gi},'_anova'])
            t = dataset2table(glme.Coefficients)
            t.fomula{1} = char(glme.Formula);
            writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
                'FileType','spreadsheet','Sheet',[varStr{q},groupStr{gi},'_Coefficients'])
        end
    end

    %%

    for gi = 1:2
        tmp_tbl = tbl(double(tbl.group) == gi,:) ;
        glme = fitglme(tmp_tbl,'beha ~ eeg*quarter+(1|subj)');
        t = dataset2table(glme.anova)
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
        t = dataset2table(glme.Coefficients);
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
    end
    %% plot

    for gi = 2%1:2
        for q = 2%1:2
            for ssn = 1:3
                for ssn1 = 1:3
                    tmp_idx = double(tbl.group)==gi & double(tbl.quarter)==q &double(tbl.ssn)==ssn &double(tbl.ssn1)==ssn1;
                    tmp_tbl = tbl(tmp_idx,:);
                    figure;
                    scatter(tmp_tbl.eeg,tmp_tbl.beha)
                    lsline;
                    [r,p] = corr(tmp_tbl.eeg,tmp_tbl.beha);
                    glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)');
                    t = dataset2table(glme.anova);

                    text(0.1,0.9,sprintf('r=%.3f,p= %.3f\nglme F=%.3f,p=%.3f',r,p,t.FStat(2),t.pValue(2)),'sc')
                    title([groupStr{gi} varStr{q}],sprintf('ssN=%d,ssN+1= %d',ssn,ssn1))

                end
            end
        end
    end


    %% separate group: frontal
    for gi = 1:2
        tmp_tbl = tbl(double(tbl.group) == gi,:) ;
        glme = fitglme(tmp_tbl,'beha ~ eeg*ss+(1|subj)');

        t = dataset2table(glme.anova);
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[groupStr{gi},'_anova'])
        t = dataset2table(glme.Coefficients)
        t.fomula{1} = char(glme.Formula);
        writetable(t,fullfile(Dir.ana,'TableOutput',['RespVariabilityGuidedPower_seq_N+1glme_SampledAT',num2str(Nsamps),'_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.xls']),...
            'FileType','spreadsheet','Sheet',[groupStr{gi},'_Coefficients'])
    end

    gi = 1;
    for s = 1:3
        tmp_tbl = tbl(double(tbl.group) == gi & double(tbl.ss) == s,:) ;
        glme = fitglme(tmp_tbl,'beha ~ eeg+(1|subj)');
        t = dataset2table(glme.anova)
    end

    %% get the correlation coefficients of each subj and hist plot

    % histogram of correlation Rho
    myColors = flipud(crameri('bam',2));% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
    grougAligStr = {'Right','Left'};
    figure('position',[100 100 600 400])
    for q = 1:2
        subplot(1,2,q);axis square;  hold all;
        for gi  = 1:2
            tmpsubs = find(subs.group == gi);
            [~,p,~,stats] = ttest(subsAll.corrR(tmpsubs,q));

            h = histogram(subsAll.corrR(tmpsubs,q));
            h.Normalization = 'probability';
            h.BinWidth = 0.14;
            h.FaceColor = myColors(gi,:);
            h.FaceAlpha = 0.9;
            plot([mean(subsAll.corrR(tmpsubs,q)) mean(subsAll.corrR(tmpsubs,q))],[0 1],'color',myColors(gi,:),'LineWidth',2,'HandleVisibility','off')
            text(2.2-1.1*gi,0.8,sprintf('Rho = %.3f\nt=%.3f\np=%.3f\nd=%.3f',mean(subsAll.corrR(tmpsubs,q)),stats.tstat,p,computeCohen_d(subsAll.corrR(tmpsubs,q),zeros(length(tmpsubs),1),'paired')),...
                'sc','color',myColors(gi,:),'HorizontalAlignment',grougAligStr{gi})
        end
        legend(groupStr,'Location','southoutside');
        ylabel('Proportion')
        xlabel('Spearman Rho value')
        title([chanROI{:}],varStr{q})
        set(gca,'ylim',[0 0.5],'xlim',[-0.8 0.5])
    end

    saveas(gca,fullfile(Dir.figs,['RespVarianceGuidedPowerRT_seq_corr_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))
end