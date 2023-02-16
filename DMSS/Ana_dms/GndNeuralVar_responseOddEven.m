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
condStr = {'ss1','ss2','ss4','collapse'};

% frontalROI = {'Fz','F1','F2'};
% frontalROI = {'Fz','F1','F2','FCz','FC1','FC2'};
frontalROI = {'CPz','CP1','CP2','Pz','Cz'};
% frontalROI = {'POz','Oz'};
occipROI = {'POz','Oz'};
freq.betaFreq = [15 25];% Hz
freq.alphaFreq = [8 13];% Hz

timeROI.all = [-0.4 0.5];% in s, for curve plotting
timeROI.Post = [0 0.5];% in s
timeROI.Pre = [-0.4 0];% in s

%%
clear subsAll
for sub_i = 1:subN
    subname = subs.name{sub_i};
    matName =fullfile(Dir.results,[subname,'_resp_ft_var_oddeven',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);

        freqID.beta = dsearchn(tfDat{1}.freq',freq.betaFreq');
        freqID.beta = freqID.beta(1):freqID.beta(2);

        freqID.alpha = dsearchn(tfDat{1}.freq',freq.alphaFreq');
        freqID.alpha = freqID.alpha(1):freqID.alpha(2);

        timeID.all = dsearchn(tfDat{1}.time',timeROI.all');
        timeID.all = timeID.all(1):timeID.all(2);

        timeID.Post = dsearchn(tfDat{1}.time',timeROI.Post');
        timeID.Post = timeID.Post(1):timeID.Post(2);

        timeID.Pre = dsearchn(tfDat{1}.time',timeROI.Pre');
        timeID.Pre = timeID.Pre(1):timeID.Pre(2);

        clear chanID
        [log1,chanID.frontal] = ismember(frontalROI,tfDat{1}.label);
        [log2,chanID.occip] = ismember(occipROI,tfDat{1}.label);

        for cond_i = 1:3
            for b = 1:size(tfDat,2)
                dat{sub_i}.betaAvgFrontal(b,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i,b}.powspctrm(chanID.frontal(log1),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));
                dat{sub_i}.betaAvgOccip(b,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i,b}.powspctrm(chanID.occip(log2),freqID.beta,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));

                dat{sub_i}.betaChans(b,cond_i,:) = squeeze(mean(mean(tfDat{cond_i,b}.powspctrm(:,freqID.beta,timeID.Post),2,"omitnan"),3,"omitnan"));


                dat{sub_i}.alphaAvgFrontal(b,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i,b}.powspctrm(chanID.frontal(log1),freqID.alpha,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));
                dat{sub_i}.alphaAvgOccip(b,cond_i) = squeeze(mean(mean(mean(tfDat{cond_i,b}.powspctrm(chanID.occip(log2),freqID.alpha,timeID.Post),1,"omitnan"),2,"omitnan"),3,"omitnan"));

                dat{sub_i}.alphaChans(b,cond_i,:) = squeeze(mean(mean(tfDat{cond_i,b}.powspctrm(:,freqID.alpha,timeID.Post),2,"omitnan"),3,"omitnan"));


                dat{sub_i}.behaACC(b,cond_i) = mean(PostTrlBeha{cond_i,b}.button_resp_corr,'omitnan');
                dat{sub_i}.behaRT(b,cond_i) = mean(PostTrlBeha{cond_i,b}.button_resp_rt,'omitnan');
            end



            %                 [subsAll.betaFrontalRTcorr(sub_i,cond_i),subsAll.betaFrontalRT_pval(sub_i,cond_i)]= corr(dat{sub_i}.betaAvgFrontal(:,cond_i),dat{sub_i}.behaRT(:,cond_i),'type','Spearman','rows','complete');
            %                 [subsAll.betaFrontalACCcorr(sub_i,cond_i),subsAll.betaFrontalACC_pval(sub_i,cond_i)] = corr(dat{sub_i}.betaAvgFrontal(:,cond_i),dat{sub_i}.behaACC(:,cond_i),'type','Spearman','rows','complete');
            %                 [subsAll.betaOccipRTcorr(sub_i,cond_i),subsAll.betaOccipRT_pval(sub_i,cond_i)] = corr(dat{sub_i}.betaAvgOccip(:,cond_i),dat{sub_i}.behaRT(:,cond_i),'type','Spearman','rows','complete');
            %                 [subsAll.betaOccipACCcorr(sub_i,cond_i),subsAll.betaOccipACC_pval(sub_i,cond_i)] = corr(dat{sub_i}.betaAvgOccip(:,cond_i),dat{sub_i}.behaACC(:,cond_i),'type','Spearman','rows','complete');
            %
            %                 subsAll.alphaFrontalRTcorr(sub_i,cond_i) = corr(dat{sub_i}.alphaAvgFrontal(:,cond_i),dat{sub_i}.behaRT(:,cond_i),'type','Spearman','rows','complete');
            %                 subsAll.alphaFrontalACCcorr(sub_i,cond_i) = corr(dat{sub_i}.alphaAvgFrontal(:,cond_i),dat{sub_i}.behaACC(:,cond_i),'type','Spearman','rows','complete');
            %                 subsAll.alphaOccipRTcorr(sub_i,cond_i) = corr(dat{sub_i}.alphaAvgOccip(:,cond_i),dat{sub_i}.behaRT(:,cond_i),'type','Spearman','rows','complete');
            %                 subsAll.alphaOccipACCcorr(sub_i,cond_i) = corr(dat{sub_i}.alphaAvgOccip(:,cond_i),dat{sub_i}.behaACC(:,cond_i),'type','Spearman','rows','complete');

        end
        subsAll.betaFrontalVar(sub_i,:,:) = dat{sub_i}.betaAvgFrontal;
        subsAll.behaACC(sub_i,:,:) = dat{sub_i}.behaACC;

        %             for c = 1:length(tfDat{1}.label)
        %                 tmp_var = squeeze(dat{sub_i}.betaChans(:,:,c));
        %                 tmp_beha =  dat{sub_i}.behaACC(:,cond_i);
        %                 subsAll.betaChansCorrACC(sub_i,c,:) = corr(tmp_var,tmp_beha,'rows','pairwise');
        %
        %                 tmp_beha =  dat{sub_i}.behaRT(:,cond_i);
        %                 subsAll.betaChansCorrRT(sub_i,c,:) = corr(tmp_var,tmp_beha,'rows','pairwise');
        %
        %                 tmp_var = squeeze(dat{sub_i}.alphaChans(:,:,c));
        %                 tmp_beha =  dat{sub_i}.behaACC(:,cond_i);
        %                 subsAll.alphaChansCorrACC(sub_i,c,:) = corr(tmp_var,tmp_beha,'rows','pairwise');
        %
        %                 tmp_beha =  dat{sub_i}.behaRT(:,cond_i);
        %                 subsAll.alphaChansCorrRT(sub_i,c,:) = corr(tmp_var,tmp_beha,'rows','pairwise');
        %
        %             end
    end
end
%     save(fullfile(Dir.results,['RespBetaVar_SampledAT',num2str(Nsamps),num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.all(1)),'~',num2str(timeROI.all(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']),'subsAll','-v7.3')
%%
figure('position',[200 200+b*100 800 1000]);
binStr = {'Odd','Even'};
for b = 1:2
    clear r p
    for gi = 1:2
        for cond_i = 1:3
            X=subsAll.betaFrontalVar(subs.group==gi,b,cond_i);
            Y=subsAll.behaACC(subs.group==gi,b,cond_i);
            [r(gi,cond_i),p(gi,cond_i)] = corr(X,Y,'type','Spearman');

            subplot(4,4,(gi-1)*4+cond_i+(b-1)*8);hold all;axis square
            plot(X,Y,'k.')
            lsline
            text(1,0.2,sprintf('Rho=%.3f\np=%.3f',r(gi,cond_i),p(gi,cond_i)),'sc','HorizontalAlignment','right')
            title(condStr(cond_i),[groupStr{gi} binStr{b}]);
            xlabel('Averaged N+1 acc');
            ylabel('Beta var of N');
            set(gca,'YLim',[0.6 1]);
        end

        cond_i = 4;
        X=mean(subsAll.betaFrontalVar(subs.group==gi,b,:),3);
        Y=mean(subsAll.behaACC(subs.group==gi,b,:),3);
        [r(gi,cond_i),p(gi,cond_i)] = corr(X,Y,'type','Spearman');

        subplot(4,4,(gi-1)*4+cond_i+(b-1)*8);hold all;axis square
        plot(X,Y,'k.')
        lsline
        text(1,0.2,sprintf('Rho=%.3f\np=%.3f',r(gi,cond_i),p(gi,cond_i)),'sc','HorizontalAlignment','right')
            title(condStr(cond_i),[groupStr{gi} binStr{b}]);
        xlabel('Averaged N+1 acc');
        ylabel('Beta var of N');
        set(gca,'YLim',[0.6 1]);
    end
end

saveas(gca,fullfile(Dir.figs,['RespBetaOddEven_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']))

