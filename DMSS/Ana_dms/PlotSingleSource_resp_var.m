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

freq.betaFreq = [15 25];% Hz
% timeROI.bins = [0 0.5];
timeROI.bins = [0.1 0.5];

%%
% [14 17 21]
for sub_i = 17%1:subN
    subname = subs.name{sub_i};
    outputFile =fullfile(Dir.results,[subname,'_resp_LCMV',txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(outputFile)
        load(outputFile);
        cfg = [];
        cfg.latency = timeROI.bins;
        cfg.frequency  = freq.betaFreq;
        tfDat.data = ft_freqdescriptives(cfg,tfDat.data);

        res = tfDat.sourcemodel;
        res.powspctrm = nan(length(res.inside),1);
        res.powspctrm(res.inside) = squeeze(mean(mean(tfDat.data.powspctrm,2),3));

        tmp = tfDat.sourcemodel;
        tmp.powspctrm = nan(length(tmp.inside),1);
        tmp.powspctrm(tmp.inside) = squeeze(mean(mean(tfDat.bl.powspctrm,2),3));

        res.powspctrm =  res.powspctrm - tmp.powspctrm;
        %         res.powspctrm =  (res.powspctrm - tmp.powspctrm)./tmp.powspctrm;

        %% Perform Statistical Analysis: https://macquarie-meg-research.github.io/MQ_MEG_Scripts/docs/mq_source_statistics_tutorial.html#:~:text=Source%20Statistics%20%28Cluster%20Permutation%20Tests%29%20in%20Fieldtrip%201,Plot%20...%205%20Export%20to%20Nifti%20Format%20
        % https://www.fieldtriptoolbox.org/example/source_statistics/

        thresh.max = max(res.powspctrm,[],'all');
        group_mean{1} = 0.4359; % mean of younger
        group_mean{2} = 0.41; % mean of older
        if thresh.max>group_mean{subs.group(sub_i)}
            thresh.ratio = 100;% top 20%
            thresh.cutoff= group_mean{subs.group(sub_i)};
        else
            thresh.ratio = 10;% top 10%
            thresh.cutoff = prctile(res.powspctrm, 100-thresh.ratio);
        end

        mriFile = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\anatomy\single_subj_T1.nii';
        mri= ft_read_mri(mriFile);

        cfg                 = [];
        cfg.parameter = 'powspctrm';
        cfg.interpmethod    = 'nearest';
        statint             = ft_sourceinterpolate(cfg, res, mri);

        %%

        cfg = [];
        cfg.method        = 'slice';
        cfg.zlim            = 'maxabs';
        cfg.funparameter    = 'powspctrm';
        cfg.funcolormap = 'hot';
        cfg.funcolorlim = [-thresh.max thresh.max];
        cfg.opacitylim = [thresh.cutoff thresh.max];
        cfg.opacitymap    = 'rampup';

        ft_sourceplot(cfg, statint);
        title(subs.groupStr{sub_i},[num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2))]);
        %     print(fullfile(Dir.figs,'IndividualResults',[subname,'Resp_source_Slice_',subs.groupStr{sub_i},num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2)),txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']),'-dpng','-r300');

        %%

        cfg                = [];
        cfg.method         = 'surface';
        cfg.funparameter   = 'powspctrm';
        cfg.maskparameter = cfg.funparameter;
        cfg.opacitymap    = 'rampup';
        cfg.projmethod     = 'nearest';
        %     cfg.surfdownsample = 20;
        cfg.camlight       = 'no';

        %         cfg.funcolorlim = [0 thresh.max*1.5];
        cfg.funcolorlim = [0 3];
        cfg.opacitylim = [thresh.cutoff thresh.max];

        cfg.surfinflated   = 'surface_white_both.mat';

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
            c.Location = 'southoutside';
%             c.Label.Rotation = -90;
%             c.Label.Position =  c.Label.Position+[1.5 0 0];
            c.Ticks = 0:3;
%             c.Label.String = sprintf('Top %d percentage',thresh.ratio);
            c.FontSize = 16;

%             title([subname,subs.groupStr{sub_i},num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2))]);
            saveas(gcf,fullfile(Dir.figs,'IndividualResults',[subname,'Resp_source_',subs.groupStr{sub_i},'_top_',num2str(thresh.ratio),'_',num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2)),'_',num2str(a),txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.png']));
            print(gcf,fullfile(Dir.figs,'IndividualResults',[subname,'Resp_source_',subs.groupStr{sub_i},'_top_',num2str(thresh.ratio),'_',num2str(timeROI.bins(1)),'~',num2str(timeROI.bins(2)),'_',num2str(a),txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.pdf']),'-dpdf','-r600');
        end
    end
end

