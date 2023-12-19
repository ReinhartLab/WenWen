clear
load('subs.mat');
txtCell = {'','';'_lap','_dephase'};
IsdePhase = 1;

%% https://mailman.science.ru.nl/pipermail/fieldtrip/2016-July/010674.html

atlasFile = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\brainnetome\BNA_MPM_thr25_1.25mm.nii';
atlas = ft_read_atlas(atlasFile);
atlas = ft_convert_units(atlas,'cm');
% find(contains(atlas.tissuelabel,'A10m'))
% atlas.tissuelabel([13 14 27 28])
% % ROI.roi = [1 2 9 10];%A8m, A6m,http://atlas.brainnetome.org/bnatlas.html
ROI.roi = [13 14];%A10m,http://atlas.brainnetome.org/bnatlas.html
% ROI.roi = [183 184];%A24cd_l,http://atlas.brainnetome.org/bnatlas.html
% ROI.roi = [11 12 13 14];%A9m A10m,http://atlas.brainnetome.org/bnatlas.html

% 
% atlasFile = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\aal\ROI_MNI_V4.nii';
% atlas = ft_read_atlas(atlasFile);
% atlas = ft_convert_units(atlas,'cm');
% ROI.roi = [4 17];

ROI.coord = [];
for i = 1:atlas.dim(1)
    for j = 1:atlas.dim(2)
        for k = 1:atlas.dim(3)
            if ismember(atlas.tissue(i,j,k),ROI.roi)
                tmp_coord = atlas.transform*[i;j;k;1];
                ROI.coord = [ROI.coord;tmp_coord'];
            end
        end
    end
end

% ROI.r = 1.5;% radius in cm
ROI.r = 1;% radius in cm

%%
for sn = 1:height(subs)
    subname = subs.name{sn};
    if subs.excluded(sn)==1
        continue
    end
%     outFile = fullfile(Dir.results,[subname,'_PLVsource_3_','0.1~0.5',txtCell{IsdePhase+1,2},'.mat']);% 3-9 Hz

    outFile =fullfile(Dir.results,[subname,'_PLVsource','0.1~0.5',txtCell{IsdePhase+1,2},'.mat']);%4-8Hz
%         outFile = fullfile(Dir.results,[subname,'_PLVsource','0.15~0.45',txtCell{IsdePhase+1,2},'.mat']);

    if isfile(outFile)
        x = load(outFile);
    end

    if sn == 1 % getting the indices once
        ROI.in = [];
        for d = 1:length(x.source_conn_full.cor.inside)

            if x.source_conn_full.cor.inside(d)
                dis = sum(power([x.source_conn_full.cor.pos(d,:)-ROI.coord(:,1:3)],2),2);% calculate distance
%                 wi = 0; % not within ROI
%                 for c = 1:size(ROI.coord,1)
%                     if  all(x.source_conn_full.cor.pos(d,:) == ROI.coord(c,1:3))
%                         wi = 1;
%                     end
%                 end

                
                if any(dis<=ROI.r^2)
                    ROI.in = [ROI.in;d]; % indice
                end
            end
        end

    end

    x.source_conn_full.cor.plvspctrm = mean(x.source_conn_full.cor.plvspctrm(ROI.in,:),'omitnan');
    x.source_conn_full.incor.plvspctrm = mean(x.source_conn_full.incor.plvspctrm(ROI.in,:),'omitnan');

    subsAll.subs{sn,1} = x.source_conn_full.cor;
    subsAll.subs{sn,2} = x.source_conn_full.incor;

    cfg           = [];
    cfg.operation = 'x1-x2';
    cfg.parameter = 'plvspctrm';
    subsAll.subs{sn,3} = ft_math(cfg, x.source_conn_full.incor, x.source_conn_full.cor);
end
%%
idx = cellfun(@isempty,subsAll.subs(:,1));
subsAll.subs = subsAll.subs(~idx,:);
subs = subs(~idx,:);

condStr = {'Correct','Incorrect','Incor-Corr'};
%%
clear grandavg
for cond_i = 1:3
    cfg = [];
    cfg.parameter          = 'plvspctrm';
    grandavg{cond_i} = ft_sourcegrandaverage(cfg, subsAll.subs{:,cond_i});

    cfg.keepindividual     = 'yes';
    tmp= ft_sourcegrandaverage(cfg, subsAll.subs{:,cond_i});
    grandavg{cond_i}.indv = tmp.plvspctrm;
end
%% Perform Statistical Analysis: https://macquarie-meg-research.github.io/MQ_MEG_Scripts/docs/mq_source_statistics_tutorial.html#:~:text=Source%20Statistics%20%28Cluster%20Permutation%20Tests%29%20in%20Fieldtrip%201,Plot%20...%205%20Export%20to%20Nifti%20Format%20
% https://www.fieldtriptoolbox.org/example/source_statistics/
threshP.ttest = 0.05;
threshP.cluster = 0.05;
cfg                     = [];
cfg.dim                 = subsAll.subs{1}.dim;
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.parameter           = 'plvspctrm';
cfg.correctm            = 'cluster';
cfg.computecritval      = 'yes';
cfg.numrandomization    = 1000;
cfg.tail                = 0;
cfg.alpha               = threshP.ttest; % alpha level of the permutation test
cfg.clusteralpha        = threshP.cluster; % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clustertail         = 0;
% Design Matrix
nsubj                   = height(subs);
cfg.design(1,:)         = [1:nsubj 1:nsubj];
cfg.design(2,:)         = [ones(1,nsubj) ones(1,nsubj)*2];
cfg.uvar                = 1;% row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar                = 2; % row of design matrix that contains independent variable (the conditions)

stat                    = ft_sourcestatistics(cfg,subsAll.subs{:,2},subsAll.subs{:,1}); % incor minus cor

stat.mask(ROI.in) = nan;% mask out seed region
stat.stat(ROI.in) = nan;% mask out seed region
stat.posclusterslabelmat(ROI.in) = nan;% mask out seed region
stat.negclusterslabelmat(ROI.in) = nan;% mask out seed region

%% Interpolate onto SPM T1 Brain

mriFile = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\anatomy\single_subj_T1.nii';
mri= ft_read_mri(mriFile);

cfg                 = [];
cfg.parameter       = 'stat';
cfg.interpmethod    = 'nearest';
statint             = ft_sourceinterpolate(cfg, stat, mri);

cfg.parameter       = 'mask';
maskint             = ft_sourceinterpolate(cfg, stat,mri);
statint.mask        = maskint.mask;

%%
cfg                 = [];
cfg.method        = 'ortho';
cfg.zlim            = 'maxabs';
cfg.funparameter    = 'stat';
cfg.maskparameter   = 'mask';
cfg.location        = 'max';
cfg.funcolormap = 'hot';
% cfg.atlas = 'D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\aal\ROI_MNI_V4.nii';
% cfg.queryrange = 3;
ft_sourceplot(cfg,statint);

% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

print(fullfile(Dir.figs,['source_int_2D_',strrep(num2str(ROI.roi), ' ', ''),num2str(x.source_conn_full.toi(1)),'~',num2str(x.source_conn_full.toi(2)),txtCell{IsdePhase+1,2},'.png']),'-dpng','-r300');

%%
cfg = [];
cfg.method        = 'slice';
cfg.zlim            = 'maxabs';
cfg.funparameter    = 'stat';
cfg.maskparameter   = 'mask';
cfg.funcolormap = 'hot';
% cfg.funcolorlim = [0 4];
ft_sourceplot(cfg, statint);
title(condStr{3});
print(fullfile(Dir.figs,['source_Slice',strrep(num2str(ROI.roi), ' ', ''),num2str(x.source_conn_full.toi(1)),'~',num2str(x.source_conn_full.toi(2)),txtCell{IsdePhase+1,2},'.png']),'-dpng','-r300');

%%

cfg                = [];
cfg.method         = 'surface';
cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.projmethod     = 'nearest';
cfg.camlight       = 'no';
cfg.surfinflated   = 'surface_white_both.mat';
    %         cfg.surfinflated   = 'surface_pial_both.mat';
% viewAng = {[0 90],[-90 0],[90 0],[180 0]};

viewAng = {[180 100],[-95 20],[100 20],[180 25]};

for a = 1:length(viewAng)

    ft_sourceplot(cfg, statint);
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
    light ('Position',[0 0 50])
    % light ('Position',[0 -50 0])

    material dull;
    drawnow;
    view(viewAng{a});
    rotate3d on
    set(gca,'FontSize',14);
    c = colorbar;
    c.Label.Rotation = -90;
    c.Label.Position = [4 0 0];
    % c.TickLabels = [-6:2:6];
    c.Label.String = sprintf('Cluster permuated p<%.3f&%.3f',threshP.ttest,threshP.cluster);
    title(condStr{3});
    print(fullfile(Dir.figs,['source_int_3D_',strrep(num2str(ROI.roi), ' ', ''),num2str(x.source_conn_full.toi(1)),'_',num2str(x.source_conn_full.toi(2)),'_',num2str(a),txtCell{IsdePhase+1,2},'.png']),'-dpng','-r300');
end


%%

% cmap = colormap("jet");
% Trange = [-6 6];
% tmp = linspace(Trange(1),Trange(2),length(cmap));
% p_alpha = 0.005;
% tvalue = norminv(1-p_alpha);
% cmap(tmp<tvalue & tmp>-tvalue,:) = repmat([1 1 1],sum(tmp<tvalue & tmp>-tvalue),1);

