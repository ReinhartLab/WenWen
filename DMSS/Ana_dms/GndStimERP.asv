clear;
load('subs.mat');
IsCorrect = 1;
ssN = [1 2 4];
for sn = 1:height(subs)
    if subs.rawEEG(sn)==1

        subname = subs.name{sn};
        set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
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
        tmptrlN = [];
        for ss = 1:3
            sEEG = pop_select(EEG,'trial',find(M.ss_num == ssN(ss)),'time',toi{ss});
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
%%
clear grandAvg
cfg = [];
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
for gi = 1:2
    for ss = 1:3
        grandAvg{gi,ss} = ft_timelockgrandaverage(cfg, subsAll{subs.rawEEG==1&subs.group==gi,ss});
        grandAvg{gi,ss}.avg = squeeze(mean(grandAvg{gi,ss}.individual));
    end
end
%%
myFigBasic

condStr = {'ss1','ss2','ss4'};
groupStr = {'Young','Old'};
chans = {'Pz','POz','Oz'};
[~,chanID] = ismember(chans,{EEG.chanlocs.labels});

figure('Position',[100 100 1200 900]);
for gi = 1:2
    for ss = 1:3
        subplot(3,4,(ss-1)*4+gi*2-1);hold all;

        dat = nan(1,length(grandAvg{1,3}.time));
        dat(1-length(grandAvg{gi,ss}.time)+end:end) = mean(grandAvg{gi,ss}.avg(chanID,:));
        plot(grandAvg{1,3}.time,dat)
        ytickformat('%,.2f')

        xlabel({'\bfTime (s)'});ylabel('\bfuV');
        set(gca,'xtick',[0 3],'xlim',toi{3},'ytick',-4:4:4,'ylim',[-4 4]);
            plot(get(gca,'xlim'),[0 0],'k','HandleVisibility','off')
            plot([0 0],get(gca,'ylim'),'k','HandleVisibility','off')
        plot([3 3],get(gca,'ylim'),'k','HandleVisibility','off')
        %     plot([3.2 3.2],get(gca,'ylim'),'k:','HandleVisibility','off')
        title([groupStr{gi} '-' condStr{ss}],[chans{:}]);
    end
end

for gi = 1:2
    for ss = 1:3
        subplot(3,4,(ss-1)*4+gi*2);hold all;axis square
        cfg = [];
        cfg.xlim = [0 3];
        cfg.zlim = [-0.4 0.4];
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
saveas(gcf,fullfile(Dir.figs,'ERP.bmp'))