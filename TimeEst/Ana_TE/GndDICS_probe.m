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

        %                         load(fullfile(Dir.results,[subname,'_dicsProbe',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']))
        load(fullfile(Dir.results,[subname,'_dicsProbe05',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']))

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
end

for gi = 1:2
    for condi = 1:3
        cfg = [];
        cfg.keepindividual     ='yes';
        gndTF{gi,condi} = ft_sourcegrandaverage(cfg, tfAll.subs{gi}{:,condi});
    end
end

% condition difference
for gi =1:2
    gndTFdiff{gi,1} = gndTF{gi,1};
    gndTFdiff{gi,1}.pow = gndTF{gi,2}.pow - gndTF{gi,1}.pow;

    gndTFdiff{gi,2} = gndTF{gi,1};
    gndTFdiff{gi,2}.pow = gndTF{gi,3}.pow - gndTF{gi,1}.pow;

    gndTFdiff{gi,3} = gndTF{gi,1};
    gndTFdiff{gi,3}.pow = gndTF{gi,3}.pow - gndTF{gi,2}.pow;
end
%%
for gi = 1:2
    for condi = 1:3
        gndTF{gi,condi}.avg.pow = mean(gndTF{gi,condi}.pow,'omitnan');
        gndTF{gi,condi}.avg.t = nan(1,length( gndTF{gi,condi}.avg.pow));
        gndTF{gi,condi}.avg.p = nan(1,length( gndTF{gi,condi}.avg.pow));

        gndTFdiff{gi,condi}.avg.pow = mean(gndTFdiff{gi,condi}.pow,'omitnan');
        gndTFdiff{gi,condi}.avg.t = nan(1,length( gndTFdiff{gi,condi}.avg.pow));
        gndTFdiff{gi,condi}.avg.p = nan(1,length( gndTFdiff{gi,condi}.avg.pow));

        for di = 1:length(gndTF{gi,condi}.avg.pow)
            if sum(isnan(gndTF{gi,condi}.pow(:,di)))==length(mySubs{gi}) %all nans
                continue
            else
                [~,gndTF{gi,condi}.avg.p(di),~,stats] =  ttest(gndTF{gi,condi}.pow(:,di));% will ignore nan
                gndTF{gi,condi}.avg.t(di) = stats.tstat;
            end
        end

        for di = 1:length(gndTFdiff{gi,condi}.avg.pow)
            if sum(isnan(gndTFdiff{gi,condi}.pow(:,di)))==length(mySubs{gi}) %all nans
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

cmap = colormap("jet");
Trange = [-6 6];
tmp = linspace(Trange(1),Trange(2),length(cmap));
p_alpha = 0.005;
tvalue = norminv(1-p_alpha);
cmap(tmp<tvalue & tmp>-tvalue,:) = repmat([1 1 1],sum(tmp<tvalue & tmp>-tvalue),1);

mri = ft_read_mri('D:\intWM-E\toolbox\fieldtrip-20211209\template\anatomy\single_subj_T1.nii');
if IsBL
    for gi = 1:2
        for condi = 1:3

            sourcePlot = gndTF{gi,condi};
            sourcePlot.avg.pow = sourcePlot.avg.t;

            sourceMask = gndTF{gi,condi};
            sourceMask.avg.pow = sourceMask.avg.pFDR;

            cfg              = [];
            cfg.parameter    = 'avg.pow';
            cfg.interpmethod = 'nearest';
            sourcePlot  = ft_sourceinterpolate(cfg, sourcePlot, mri);% variable to plot
            sourceMask  = ft_sourceinterpolate(cfg, sourceMask, mri);% mask based on pFDR
            sourcePlot.mask = sourceMask.pow < p_alpha;

%             cfg               = [];
%             cfg.method        = 'ortho';
%             cfg.funparameter  = 'pow';
%             cfg.funcolormap = 'hot';
%             cfg.funcolorlim = [0 6];
%             cfg.maskparameter = 'mask';
%             cfg.atlas = ft_read_atlas('D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\aal\ROI_MNI_V4.nii');
%             sourcePlot.coordsys = 'mni';
%             cfg.location = [-18, -36, 2];% left hippo
%             cfg.queryrange = 3;
%             source_atlas = ft_sourceplot(cfg, sourcePlot);

%                             cfg = [];
%                             cfg.method        = 'slice';
%                             cfg.funparameter  = 'pow';
%                             cfg.funcolormap = 'hot';
%                             cfg.funcolorlim = [-6 6];
%                             cfg.opacitymap    = 'vdown';
%                         cfg.opacitylim = [-0.5 0.5];
% 
%                                 cfg.maskparameter = 'mask';
%                             ft_sourceplot(cfg, sourcePlot);

            viewPara = {[0 90]};
% %                     viewPara = {[0 90],[-90 0],[90 0]};

cfg               = [];
cfg.method        = 'surface';
cfg.maskparameter = 'mask';
cfg.funparameter  = 'pow';
cfg.funcolormap = cmap;
cfg.funcolorlim = Trange;
cfg.opacitylim = [0 0.2];
cfg.opacitymap    = 'rampup';
            cfg.surfdownsample = 50;
            sourcePlot.coordsys = 'mni';
            for a = 1:length(viewPara)
                ft_sourceplot(cfg,sourcePlot);
                view(viewPara{a}(1),viewPara{a}(2))
                rotate3d on
                c = colorbar;
                c.TickLabels = [-6:2:6];
                c.Label.String = sprintf('FDR corrected t-value(%.3f)',p_alpha);
                title(condStr{condi},groupStr{gi});
                %                 saveas(gcf,fullfile('Figs',['ProbeDICS0106',condStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
%                 saveas(gcf,fullfile('Figs',['ProbeDICS',condStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
            end
        end
    end
end
%% plot condition difference
cmap = colormap("jet");
Trange = [-6 6];
tmp = linspace(Trange(1),Trange(2),length(cmap));
p_alpha = 0.005;
tvalue = norminv(1-p_alpha);
cmap(tmp<tvalue & tmp>-tvalue,:) = repmat([1 1 1],sum(tmp<tvalue & tmp>-tvalue),1);

for gi = 1:2
    for condi = 1:3

        sourcePlot = gndTFdiff{gi,condi};
        sourcePlot.avg.pow = sourcePlot.avg.t;

        sourceMask = gndTFdiff{gi,condi};
        sourceMask.avg.pow = sourceMask.avg.pFDR;

        %----to visualize mask
        %         tmm = ones(1,4050);
        %                 tmm(tpid' & idx) =0;% left temporal cortex mask
        %                 tmm(pfcid' & idx) =0;% left frontal cortext mask
%         sourceMask.avg.pow = ~(PFC.inside | tp.inside);


        cfg              = [];
        cfg.parameter    = 'avg.pow';
        cfg.interpmethod = 'nearest';
        sourcePlot  = ft_sourceinterpolate(cfg, sourcePlot, mri);% variable to plot
        sourceMask  = ft_sourceinterpolate(cfg, sourceMask, mri);% mask based on pFDR
        sourcePlot.mask = sourceMask.pow < p_alpha ;

        %                     cfg               = [];
        %                     cfg.method        = 'ortho';
        %                     cfg.funparameter  = 'pow';
        %                     cfg.funcolormap = 'hot';
        %                     cfg.maskparameter = 'mask';
        %         cfg.opacitylim = [-0.5 0.5];
        %         cfg.opacitymap    = 'vdown';
        %                     cfg.funcolorlim = [0 4];
        %                     cfg.atlas = ft_read_atlas('D:\intWM-E\toolbox\fieldtrip-20211209\template\atlas\aal\ROI_MNI_V4.nii');
        %                     sourcePlot.coordsys = 'mni';
        % %                     cfg.location = [-42 24 26];% left IFG MNI coordinate,https://academic.oup.com/cercor/article/20/4/859/305717?login=true
        %         %             cfg.location = [47 34 27];% righ MFG MNI coordinate,https://academic.oup.com/cercor/article/20/4/859/305717?login=true
        %                cfg.location = [-18, -36, 2];% left hippo
        %         cfg.queryrange = 3;
        %                     source_atlas = ft_sourceplot(cfg, sourcePlot);
        %                     title(condDiffStr{condi})

        %             cfg = [];
        %             cfg.method        = 'slice';
        %             cfg.funparameter  = 'pow';
        %             cfg.funcolormap = 'hot';
        % cfg.opacitylim = [-0.5 0.5];
        %         cfg.opacitymap    = 'vdown';
        %             cfg.maskparameter = 'mask';
        %             cfg.funcolorlim = [0 4];
        %             ft_sourceplot(cfg, sourcePlot);
        %           title(condDiffStr{condi});

%         viewPara = {[0 90],[-90 0],[90 0],[180 0]};
                        viewPara = {[180 0]};

        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = 'mask';
        cfg.funcolormap = cmap;
        cfg.funcolorlim = Trange;
        cfg.opacitylim = [0 0.2];
        cfg.opacitymap    = 'rampup';
        cfg.surfdownsample = 50;

        for a = 1:length(viewPara)
            ft_sourceplot(cfg,sourcePlot);
            % [caz,cel] = view
            view(viewPara{a})
            c = colorbar;
            c.TickLabels = [-6:2:6];
            c.Label.String = sprintf('FDR corrected t-value(%.3f)',p_alpha);
            rotate3d on
            title(condDiffStr{condi},groupStr{gi});
%             saveas(gcf,fullfile('Figs',['ProbeDICS' ,condDiffStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
            %                    saveas(gcf,fullfile('Figs',['ProbeDICS0106' ,condDiffStr{condi},groupStr{gi},num2str(a),txtCell{IsBL+1,4},txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.bmp']))
        end
    end
end
%% creating seeds: based on young group, attend minus ignore,find the max theta diff at left temporal and frontal cortex
clear PFC tp
idx = gndTFdiff{1,3}.avg.pFDR<p_alpha;
dPos = gndTFdiff{1,3}.pos;

tpid = dPos(:,1)<0 &  dPos(:,2)<2 & dPos(:,3)>-3 & dPos(:,3)<2;
[~,tmp] = max(gndTFdiff{1,3}.avg.t(tpid' & idx));
pl =find(tpid' & idx);
tp.roi = dPos(pl(tmp),:)*10

pfcid = dPos(:,1)<0 &  dPos(:,2)>-3 & dPos(:,3)>0;
[~,tmp] = max(gndTFdiff{1,3}.avg.t(pfcid' & idx));
pl =find(pfcid' & idx);
PFC.roi = dPos(pl(tmp),:)*10


sphereR = 35;% in mm

PFC.inside = zeros(size(dPos,1),1);
tp.inside = zeros(size(dPos,1),1);

for t = 1:size(dPos,1)
    if sum((dPos(t,:)*10-PFC.roi).^2) <=sphereR^2
        PFC.inside(t) = 1;
    end
    if sum((dPos(t,:)*10-tp.roi).^2) <=sphereR^2
        tp.inside(t) = 1;
    end

end
PFC.inside = PFC.inside & idx';
tp.inside = tp.inside & idx';

C = intersect(find(PFC.inside),find(tp.inside));
PFC.inside(C) = 0;
tp.inside(C) = 0;

sum(PFC.inside)
sum(tp.inside)
save('PFCTP.mat','PFC','tp')

