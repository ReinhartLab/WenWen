clear;
load('subs.mat');
subs(subs.rawEEG==0 | subs.exclude==1,:) = [];

ssN = [1 2 4];

for IsCorrect = [1 0]
    for IsBL2preDelay = [0 1]

        for sn = 1:height(subs)
            subname = subs.name{sn};
            set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
            if isfile(set_name)
                EEG = pop_loadset('filename',set_name);

                %% load behavior

                csvFile = fullfile(Dir.beha,subs.csvFile{sn});
                M = readtable(csvFile);
                M = M(:,["block_num","ss_num","type","button_resp_rt","button_resp_corr"]);
                M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});

                M(isnan(M.ss_num),:) = [];
                if sum(isnan(M.button_resp_rt)) >0
                    M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
                end

                tmpEEG = pop_select(EEG,'time',[0 2]);%to avoid multiple trial marker within long epochs
                tmp = find(contains({tmpEEG.event.type},'T'));
                if length(tmp) ==EEG.trials
                    cleanTrials = {tmpEEG.event(tmp).type};
                    cleanTrials = cellfun(@(x)(str2num(x(2:end))),cleanTrials,'UniformOutput',false);
                    cleanTrials = cell2mat(cleanTrials);
                end
                M = M(cleanTrials,:);clear tmpEEG

                %%
                if IsCorrect
                    EEG = pop_select(EEG,'trial',find(M.button_resp_corr==1));
                    M(M.button_resp_corr~=1,:) = [];
                end

                %%
                if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
                    stim1 = {'S 11','S 21', 'S 31'};
                else
                    stim1 = {'11','21','31'};
                end

                toi = {[-1.4 3.5],[-2.6 3.5],[-5 3.5]};% time-locked2delay
                bsWindow = {[-1.4 -1.2],[-2.6 -2.4],[-5 -4.8]};
                tmptrlN = [];
                for ss = 1:3
                    sEEG = pop_select(EEG,'trial',find(M.ss_num == ssN(ss)),'time',toi{ss});
                    if IsBL2preDelay
                        sEEG = pop_rmbase(sEEG,[-200 0]);% pre delay baseline
                    else
                        sEEG = pop_rmbase(sEEG,bsWindow{ss}*1000);% pre trial baseline
                    end
                    tmptrlN(ss) =sEEG.trials;

                    eeg = eeglab2fieldtrip(sEEG,'preprocessing','none');
                    cfg = [];
                    avgTrl = ft_timelockanalysis(cfg,eeg);
                    subsAll{sn,ss} = avgTrl;
                end
                if EEG.trials-sum(tmptrlN)>10
                    error
                end
            end
        end
        txtCell = {'','';'_corrTrials','_bl2preDelay'};

        save(['ERP',txtCell{IsCorrect+1,1},txtCell{IsBL2preDelay+1,2},'.mat'],'subsAll','-v7.3')
        %%
        load(['ERP',txtCell{IsCorrect+1,1},txtCell{IsBL2preDelay+1,2},'.mat'])

        empty_idx = cellfun(@isempty,subsAll(:,1));
        subsAll(empty_idx,:) = [];
        subs(empty_idx,:) = [];
        %%
        clear grandAvg

        for gi = 1:2
            for ss = 1:3
                cfg = [];
                cfg.parameter = 'avg';
                cfg.keepindividual = 'yes';
                tmp = ft_timelockgrandaverage(cfg, subsAll{subs.rawEEG==1&subs.group==gi,ss});

                cfg = [];
                cfg.parameter = 'avg';
                grandAvg{gi,ss} = ft_timelockgrandaverage(cfg, subsAll{subs.rawEEG==1&subs.group==gi,ss});
                grandAvg{gi,ss}.indv = tmp.individual;
            end
        end

        cfg = [];
        cfg.parameter = 'avg';
        cfg.operation = 'subtract';
        for ss = 1:3
            groupDif{ss} = ft_math(cfg, grandAvg{2,ss},grandAvg{1,ss});% old minus young
        end
        %%
        myFigBasic

        condStr = {'ss1','ss2','ss4'};
        groupStr = {'Young','Old','Old-Young'};
        chans = {'Pz','POz','Oz'};
        [~,chanID] = ismember(chans,{EEG.chanlocs.labels});

        figure('Position',[100 100 1300 900]);
        for gi = 1:2
            for ss = 1:3
                subplot(3,5,(ss-1)*5+gi*2-1);hold all;

                dat = nan(1,length(grandAvg{1,3}.time));
                dat(1-length(grandAvg{gi,ss}.time)+end:end) = mean(grandAvg{gi,ss}.avg(chanID,:));
                plot(grandAvg{1,3}.time,dat)
                ytickformat('%.1f')

                xlabel({'Time(0s=maintenance)'});
                set(gca,'xtick',[-4.8 -2.4 -1 0 3],'xlim',toi{3},'ytick',-2:2:6,'ylim',[-3 6]);
                plot(get(gca,'xlim'),[0 0],'k','HandleVisibility','off')
                plot([-1 -1],get(gca,'ylim'),'k--','HandleVisibility','off')
                plot([0 0],get(gca,'ylim'),'k','HandleVisibility','off')
                plot([3 3],get(gca,'ylim'),'k','HandleVisibility','off')
                %     plot([3.2 3.2],get(gca,'ylim'),'k:','HandleVisibility','off')
                title([groupStr{gi} '-' condStr{ss}],[chans{:}]);
            end
        end

        for gi = 1:2
            for ss = 1:3
                subplot(3,5,(ss-1)*5+gi*2);hold all;axis square
                cfg = [];
                cfg.xlim = [0 3];
                cfg.zlim = [-1 1];
                cfg.colormap = jet;
                cfg.highlight          =  'on';
                cfg.highlightsymbol    = '.';
                cfg.highlightsize      = 24;
                cfg.highlightchannel = chans;
                cfg.comment  = 'auto';
                cfg.markersize = 3;
                cfg.marker = 'off';
                cfg.figure = 'gca';

                cfg.layout = 'easycapM11.mat';
                ft_topoplotER(cfg,grandAvg{gi,ss});

                h = colorbar;
                % set(h,'ytick',-5:1:10);
                %         h.Label.String = '\bf\muV';
                title(groupStr{gi}, condStr{ss});
            end
        end


        for ss = 1:3
            subplot(3,5,ss*5);hold all;axis square
            cfg = [];
            cfg.xlim = [0 3];
            cfg.zlim = [-1 1];
            cfg.colormap = jet;

            cfg.comment  = 'auto';
            cfg.markersize = 3;
            cfg.marker = 'off';
            cfg.figure = 'gca';

            cfg.layout = 'easycapM11.mat';
            ft_topoplotER(cfg,groupDif{ss});

            h = colorbar;
            % set(h,'ytick',-5:1:10);
            h.Label.String = 'Voltage(\muV)';
            h.Label.Position = [5.1 0 0];
            h.Label.Rotation = -90;
            title(groupStr{3}, condStr{ss});
        end
        saveas(gcf,fullfile(Dir.figs,['ERP',txtCell{IsCorrect+1,1},txtCell{IsBL2preDelay+1,2},'.bmp']))

        %%
        cfg = [];
        cfg.layout = 'easycapM11.mat';
        cfg.interactive = 'yes';
        cfg.showoutline = 'yes';
        ft_multiplotER(cfg, grandAvg{gi,1}, grandAvg{gi,2},grandAvg{gi,3})
    end
end