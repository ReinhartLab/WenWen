clear;rng('shuffle');
load('subs.mat');

subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic
addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;

IsCorretTrials = 1;
IsBL2preDelay = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Young-Old'};
condStr = {'Load 1','Load 2','Load 4'};

freq.betaFreq = [15 25];% Hz
% timeROI.bins = [0 0.5];
timeROI.bins = [0.1 0.5];

%%
subsAll = cell(subN,1);
blAll = subsAll;
for sub_i = 1:subN
    subname = subs.name{sub_i};
    outputFile =fullfile(Dir.results,[subname,'_resp_LCMV',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(outputFile)
        load(outputFile);
        cfg = [];
        cfg.latency = timeROI.bins;
        cfg.frequency  = freq.betaFreq;
        tfDat.data = ft_freqdescriptives(cfg,tfDat.data);

        tmp = tfDat.sourcemodel;
        tmp.powspctrm = nan(length(tmp.inside),1);
        tmp.powspctrm(tmp.inside) = squeeze(mean(mean(tfDat.data.powspctrm,2),3));
        subsAll{sub_i} = tmp;

        tmp = tfDat.sourcemodel;
        tmp.powspctrm = nan(length(tmp.inside),1);
        tmp.powspctrm(tmp.inside) = squeeze(mean(mean(tfDat.bl.powspctrm,2),3));
        blAll{sub_i} = tmp;
    end
end

empty_idx = cellfun(@isempty,subsAll);
subsAll(empty_idx) = [];
blAll(empty_idx) = [];
subs(empty_idx,:) = [];

%% Perform Statistical Analysis: https://macquarie-meg-research.github.io/MQ_MEG_Scripts/docs/mq_source_statistics_tutorial.html#:~:text=Source%20Statistics%20%28Cluster%20Permutation%20Tests%29%20in%20Fieldtrip%201,Plot%20...%205%20Export%20to%20Nifti%20Format%20
% https://www.fieldtriptoolbox.org/example/source_statistics/
threshP.ttest = 0.005;
threshP.cluster = 0.005;

for gi = 1:2 % separate two groups

    cfg                     = [];
    cfg.dim                 = subsAll{1}.dim;
    cfg.method              = 'montecarlo';
    cfg.parameter           = 'powspctrm';
    cfg.correctm            = 'cluster';
    cfg.numrandomization    = 10000;
    cfg.alpha               = threshP.ttest; % alpha level of the permutation test
    cfg.clusteralpha        = threshP.cluster; % alpha level of the sample-specific test statistic that will be used for thresholding
    cfg.computecritval      = 'yes';
    cfg.statistic           = 'ft_statfun_depsamplesT';

    % Design Matrix
    nsubj                   = sum(subs.group==gi);
    cfg.design(1,:)         = [1:nsubj 1:nsubj];
    cfg.design(2,:)         = [ones(1,nsubj) ones(1,nsubj)*2];
    cfg.uvar                = 1;% row of design matrix that contains unit variable (in this case: subjects)
    cfg.ivar                = 2; % row of design matrix that contains independent variable (the conditions)

    stat                    = ft_sourcestatistics(cfg,subsAll{subs.group==gi},blAll{subs.group==gi}); % don't forget the {:}!
    %% Interpolate onto SPM T1 Brain

    mriFile = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\anatomy\single_subj_T1.nii';
    mri= ft_read_mri(mriFile);

    subsAll{1}.stat = stat.stat;
    subsAll{1}.mask = stat.mask;

    cfg                 = [];
    cfg.parameter = 'stat';
    cfg.interpmethod    = 'nearest';
    %         statint             = ft_sourceinterpolate(cfg, stat, mri);
    statint             = ft_sourceinterpolate(cfg, subsAll{1}, mri);

    cfg.parameter       = 'mask';
    maskint             = ft_sourceinterpolate(cfg, subsAll{1},mri);
    statint.mask        = maskint.mask;

    %%
    threshP.tvalue = norminv(1-threshP.ttest);

    cfg = [];
    cfg.method        = 'slice';
    cfg.zlim            = 'maxabs';
    cfg.funparameter    = 'stat';
    cfg.maskparameter   = 'mask';
    cfg.funcolormap = 'hot';
    %     cfg.funcolorlim = [threshP.tvalue 6];
    cfg.funcolorlim = [0 8];
    ft_sourceplot(cfg, statint);
    title(groupStr{gi},[num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2))]);
    print(fullfile(Dir.figs,['Resp_source_Slice_',groupStr{gi},num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2)),txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']),'-dpng','-r450');

    %%

    cfg                = [];
    cfg.method         = 'surface';
    cfg.maskparameter  = 'mask';
    cfg.funparameter   = 'stat';
    cfg.projmethod     = 'nearest';
    %     cfg.surfdownsample = 20;
    cfg.camlight       = 'no';
    cfg.colormap = 'hot';
    %     cfg.funcolorlim = [threshP.tvalue 6];
    cfg.funcolorlim = [0 8];
    cfg.surfinflated   = 'surface_white_both.mat';
    %         cfg.surfinflated   = 'surface_pial_both.mat';

    % viewAng = {[0 90],[-90 0],[90 0],[180 0]};

    viewAng = {[180 100],[-95 20],[100 20],[180 25]};

    for a = 1%:length(viewAng)

        ft_sourceplot(cfg, statint);
        %         ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        %         colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        light ('Position',[0 0 50])
        %         light ('Position',[0 -50 0])

        material dull;
        drawnow;
        view(viewAng{a});
        rotate3d on
        set(gca,'FontSize',14);
        c = colorbar;
        c.Label.Rotation = -90;
        c.Label.Position =  c.Label.Position+[1.5 0 0];
        % c.TickLabels = [-6:2:6];
        c.Label.String = 'T-value of cluster based permutation';
        c.FontSize = 16;

        %         c.Label.String = sprintf('Cluster permuated p<%.3f&%.3f',threshP.ttest,threshP.cluster);
%         title([groupStr{gi},num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2))]);
        saveas(gcf,fullfile(Dir.figs,['Resp_source_int_3D_',groupStr{gi},num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2)),'_',num2str(a),txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']));
        print(gcf,fullfile(Dir.figs,['Resp_source_int_3D_',groupStr{gi},num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2)),'_',num2str(a),txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r600');
    end
end

%  cmap = colormap("jet");
% Trange = [-6 6];
% tmp = linspace(Trange(1),Trange(2),length(cmap));
% p_alpha = 0.005;
% tvalue = norminv(1-p_alpha);
% cmap(tmp<tvalue & tmp>-tvalue,:) = repmat([1 1 1],sum(tmp<tvalue & tmp>-tvalue),1);

