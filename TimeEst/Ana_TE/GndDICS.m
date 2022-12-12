clear;rng('shuffle');
load('subs.mat');

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;

IsBL = 1;% whether to apply baseline

IsLap = 1;
IsdePhase = 1;
IsCorretTrials = 1;
txtCell = {'','','','';'_lap','_dephase','_corrTrials','_BL'};
%% participants

tmp = find(subs.group==1 & ~ismember(subs.name,{'S05','S28','S27','SE08','SE16'}));

% S05: <20% EEG trials remained with a threshold of 75uV
% S27: quit after 4 blocks
% S28: bad data
mySubs{1} = tmp;% YOUNG

tmp = find(subs.group==2 & ~ismember(subs.name,{'S05','S28','S27','SE01','SE05','SE08','SE16'}));
mySubs{2} = tmp;% OLD

%%
clear tfAll gndTF
for gi = 1:2
    for sub_i = 1:length(mySubs{gi})
        sn = mySubs{gi}(sub_i);
        subname = subs.name{sn};

        load(fullfile(Dir.results,[subname,'_dics',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']))

        if IsBL
            cfg = [];
            %                 cfg.operation='subtract';
            cfg.operation = '(x1 ./ x2) - 1';
            cfg.parameter = 'avg.pow';
            for condi = 1:3
                source.DataBL{condi} = ft_math(cfg,source.Data{condi},source.BL{condi});
            end
            tfAll.subs{gi}(sub_i,:) = source.DataBL;% baseline correction
        else
            tfAll.subs{gi}(sub_i,:) = source.Data;% using unbaseline corrected data
        end
    end

    for condi = 1:3
        cfg = [];
        cfg.keepindividual     ='yes';
        gndTF{gi,condi} = ft_sourcegrandaverage(cfg, tfAll.subs{gi}{:,condi});
    end
end

% condition difference
for gi =1:2
    gndTFdiff{gi,1} = gndTF{1};
    gndTFdiff{gi,1}.pow = gndTF{gi,2}.pow - gndTF{gi,1}.pow;

    gndTFdiff{gi,2} = gndTF{1};
    gndTFdiff{gi,2}.pow = gndTF{gi,3}.pow - gndTF{gi,1}.pow;

    gndTFdiff{gi,3} = gndTF{1};
    gndTFdiff{gi,3}.pow = gndTF{gi,3}.pow - gndTF{gi,2}.pow;
end
%%
for gi = 1:2
    for condi = 1:3
        gndTF{gi,condi}.avg.pow = mean(gndTF{gi,condi}.pow,'omitnan');
        gndTF{gi,condi}.avg.t = nan(1,length( gndTF{gi,condi}.avg.pow));
        gndTF{gi,condi}.avg.p = gndTF{gi,condi}.avg.t;

        gndTFdiff{gi,condi}.avg.pow = mean(gndTFdiff{gi,condi}.pow,'omitnan');
        gndTFdiff{gi,condi}.avg.t = nan(1,length(gndTFdiff{gi,condi}.avg.pow));
        gndTFdiff{gi,condi}.avg.p = gndTFdiff{gi,condi}.avg.t;

        for di = 1:length(gndTF{gi,condi}.avg.pow)
            if sum(isnan(gndTF{gi,condi}.pow(:,di)))>=length(mySubs{gi}) %all nans
                continue
            else
                [~,gndTF{gi,condi}.avg.p(di),~,stats] =  ttest(gndTF{gi,condi}.pow(:,di));% will ignore nan
                gndTF{gi,condi}.avg.t(di) = stats.tstat;
            end
        end

        for di = 1:length(gndTFdiff{gi,condi}.avg.pow)
            if sum(isnan(gndTFdiff{gi,condi}.pow(:,di)))>=length(mySubs{gi}) %all nans
                continue
            else
                [~,gndTFdiff{gi,condi}.avg.p(di),~,stats] =  ttest(gndTFdiff{gi,condi}.pow(:,di));% will ignore nan
                gndTFdiff{gi,condi}.avg.t(di) = stats.tstat;
            end
        end

        if IsBL
            gndTF{gi,condi}.avg.pFDR = mafdr(gndTF{gi,condi}.avg.p);
        end

        gndTFdiff{gi,condi}.avg.pFDR = mafdr(gndTFdiff{gi,condi}.avg.p);
    end
end

%%

condStr = {'Fixation','Ignore','Attend'};
condDiffStr = {'Ignore-Fixation','Attend-Fixation','Attend-Ignore'};
groupStr = {'Young','Older'};

mri = ft_read_mri('D:\intWM-E\toolbox\fieldtrip-20211209\template\anatomy\single_subj_T1.nii');

cmap = colormap("jet");
Trange = [-6 6];
tmp = linspace(Trange(1),Trange(2),length(cmap));
p_alpha = 0.05;
tvalue = norminv(1-p_alpha);
cmap(tmp<tvalue & tmp>-tvalue,:) = repmat([1 1 1],sum(tmp<tvalue & tmp>-tvalue),1);

for gi = 1:2
    for condi = 1:3

        sourcePlot = gndTF{gi,condi};
        sourcePlot.avg.pow = sourcePlot.avg.t;

        sourceMask = gndTF{gi,condi};
        if IsBL
%             sourceMask.avg.pow = sourceMask.avg.pFDR;
            sourceMask.avg.pow = sourceMask.avg.p;
        end

        cfg              = [];
        cfg.parameter    = 'avg.pow';
        cfg.interpmethod = 'nearest';
        sourcePlot  = ft_sourceinterpolate(cfg, sourcePlot, mri);% variable to plot
        sourceMask  = ft_sourceinterpolate(cfg, sourceMask, mri);% mask based on pFDR
%         sourcePlot.mask = sourceMask.pow < p_alpha;
        sourcePlot.mask = sourceMask.pow < p_alpha & sourcePlot.pow > 0;% postive activation only

        %             cfg               = [];
        %             cfg.method        = 'ortho';
        %             cfg.funparameter  = 'pow';
        %             cfg.funcolormap = 'hot';
        %             cfg.maskparameter = 'mask';
        %             cfg.atlas = ft_read_atlas('D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\aal\ROI_MNI_V4.nii');
        %             sourcePlot.coordsys = 'mni';
        %             cfg.location = [-42 24 26];% left IFG MNI coordinate,https://academic.oup.com/cercor/article/20/4/859/305717?login=true
        %             %             cfg.location = [47 34 27];% righ MFG MNI coordinate,https://academic.oup.com/cercor/article/20/4/859/305717?login=true
        %             cfg.queryrange = 3;
        %             source_atlas = ft_sourceplot(cfg, sourcePlot);
        %
%                     cfg = [];
%                     cfg.method        = 'slice';
%                     cfg.funparameter  = 'pow';
%                     cfg.funcolorlim = [-6 6];
%                     cfg.funcolormap = 'hot';
%                     cfg.maskparameter = 'mask';
%                     cfg.opacitylim    = [-0.1 0.1];
%                     cfg.opacitymap    = 'vdown';
%                     ft_sourceplot(cfg, sourcePlot);
%                     title([groupStr{gi},condStr{condi}]);
%                     c = colorbar;
%         %             c.TickLabels = [0:5:30];
%                     c.Label.String = sprintf('uncorrected t-value (p<%.3f)',p_alpha);
        %             saveas(gcf,fullfile('Figs',['DotsDICS_slice' ,condStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))

        viewPara = {[0 90]};

        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.funcolormap = cmap;
        cfg.funcolorlim = Trange;
        cfg.maskparameter = 'mask';
        cfg.surfdownsample = 50;
        cfg.opacitymap    = 'rampup';
        cfg.opacitylim = [0 0.1];% mask

        for a = 1:length(viewPara)
            ft_sourceplot(cfg,sourcePlot);
            view(viewPara{a}(1),viewPara{a}(2))
            rotate3d on
            c = colorbar;
            c.TickLabels = [-6:2:6];
            c.Label.String = sprintf('uncorrected t-value(p<%.3f)',p_alpha);
            title(condStr{condi},groupStr{gi});
        end
        saveas(gcf,fullfile('Figs',['DotsDICS' ,condStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
    end
end

%% plot differences between conditions
cmap = colormap("jet");
Trange = [-6 6];
tmp = linspace(Trange(1),Trange(2),length(cmap));
p_alpha = 0.05;
tvalue = norminv(1-p_alpha);
cmap(tmp<tvalue & tmp>-tvalue,:) = repmat([1 1 1],sum(tmp<tvalue & tmp>-tvalue),1);

for gi = 1:2
    for condi = 1:3
        sourcePlot = gndTFdiff{gi,condi};
        sourcePlot.avg.pow = sourcePlot.avg.t;

        sourceMask = gndTFdiff{gi,condi};
        sourceMask.avg.pow = sourceMask.avg.pFDR;
        %         sourceMask.avg.pow = sourceMask.avg.p;

        cfg              = [];
        cfg.parameter    = 'avg.pow';
        cfg.interpmethod = 'nearest';
        sourcePlot  = ft_sourceinterpolate(cfg, sourcePlot, mri);% variable to plot
        sourceMask  = ft_sourceinterpolate(cfg, sourceMask, mri);% mask based on pFDR
        sourcePlot.mask = sourceMask.pow < p_alpha;

%         cfg               = [];
%         cfg.method        = 'ortho';
%         cfg.funparameter  = 'pow';
%         cfg.funcolormap = 'hot';
%         %          cfg.colorbar      = 'no';
%         cfg.maskparameter = 'mask';
%         cfg.opacitymap    = 'vdown';
%         cfg.opacitylim = [-0.2 0.2];% mask
%         cfg.funcolorlim = [-6 6];
%         cfg.atlas = ft_read_atlas('D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\aal\ROI_MNI_V4.nii');
%         sourcePlot.coordsys = 'mni';
%         cfg.location = [-42 24 26];% left IFG MNI coordinate,https://academic.oup.com/cercor/article/20/4/859/305717?login=true
%         %             cfg.location = [47 34 27];% righ MFG MNI coordinate,https://academic.oup.com/cercor/article/20/4/859/305717?login=true
%         cfg.queryrange = 3;
%         source_atlas = ft_sourceplot(cfg, sourcePlot);
%         %                     c = colorbar;
%         %             c.Label.String = sprintf('FDR corrected t-value(%.3f)',p_alpha);
%         title(condDiffStr{condi})

% 
%         cfg = [];
%         cfg.method        = 'slice';
%         cfg.funparameter  = 'pow';
%         cfg.funcolormap = 'hot';
%         cfg.maskparameter = 'mask';
%         cfg.opacitymap    = 'vdown';
%         cfg.opacitylim = [-0.1 0.1];% mask
%         cfg.funcolorlim = [-6 6];
%         ft_sourceplot(cfg, sourcePlot);
%         c = colorbar;
%         c.Label.String = sprintf('FDR corrected t-value(%.3f)',p_alpha);
%         %         c.Ticks = [0:1:6];
%         title(condDiffStr{condi});
%         saveas(gcf,fullfile('Figs',['DotsDICS_slice' ,condDiffStr{condi},groupStr{gi},txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))


        %         viewPara = {[0 90],[-90 0],[90 0],[180 0]};
        viewPara = {[0 90]};

        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = 'mask';
        cfg.funcolormap = cmap;
        cfg.funcolorlim = Trange;
        cfg.opacitylim = [-0.2 0.2];
        cfg.opacitymap    = 'vdown';
        cfg.surfdownsample = 50;

        for a = 1:length(viewPara)
            ft_sourceplot(cfg,sourcePlot);
            %             [caz,cel] = view
            view(viewPara{a}(1),viewPara{a}(2))
            rotate3d on
            c = colorbar;
            c.TickLabels = [-6:2:6];
            c.Label.String = sprintf('FDR corrected t-value(%.3f)',p_alpha);
            title(condDiffStr{condi},groupStr{gi});
        end
        saveas(gcf,fullfile('Figs',['DotsDICS' ,condDiffStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
    end
end

%% Perform Statistical Analysis

% cfg                     = [];
% cfg.dim                 = gndTF{1}.dim;
% cfg.method              = 'montecarlo';
% cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.parameter           = 'pow';
% cfg.correctm            = 'cluster';
% cfg.computecritval      = 'yes';
% cfg.numrandomization    = 4000;
% cfg.tail                = 0;
%
% % Design Matrix
% nsubj                   = numel(source_post_all);
% cfg.design(1,:)         = [1:nsubj 1:nsubj];
% cfg.design(2,:)         = [ones(1,nsubj) ones(1,nsubj)*2];
% % row of design matrix that contains unit variable (in this case: subjects)
% cfg.uvar                = 1;
% % row of design matrix that contains independent variable (the conditions)
% cfg.ivar                = 2;
%
% [stat]                  = ft_sourcestatistics(cfg,source_post_all{:},source_pre_all{:});

%
%
% % Show raw source level statistics (2D plot)
%
% cfg                     = [];
%
% cfg.method              = 'ortho';
%
% cfg.funparameter        = 'stat';
%
% cfg.location            = 'max';
%
% %cfg.maskparameter      = 'mask'; %turn on to show mask
%
% ft_sourceplot(cfg,stat);
%
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%
% colormap(flipud(brewermap(64,'RdBu'))) % change the colormap


