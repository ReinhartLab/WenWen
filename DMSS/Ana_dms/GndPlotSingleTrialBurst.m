clear;
load('subs.mat');
subs(subs.rawEEG==0 | subs.exclude==1,:) = [];
subN = height(subs);

myFigBasic

IsCorretTrials = 1;
IsBL2preDelay = 0;
IsLap = 1;
IsdePhase=1;

condStr = {'ss1','ss2','ss4'};

chanROI = {'Fz','F1','F2','FCz','FC1','FC2'};
SSepochwindow = {[-1.6,4.5],[-2.8,4.5],[-5.2,4.5]};
clear PlotSubj
PlotSubj.name = subs.name;
% PlotSubj.name = {'dmss06'};

[~,PlotSubj.ID] = ismember(PlotSubj.name,subs.name);
PlotSubj.group = subs.groupStr(PlotSubj.ID);

%%

w = [1 1 1];
b = [0,0,1];
k = [0,0,0];
r = [1,0,0];
bkr = createcolormap(w,b,k,r,w); % 256x3 array

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};
for sn= 1:length(PlotSubj.ID)

    subname = PlotSubj.name{sn};

    matName = fullfile(Dir.results,[subname,'_delay_ft',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);
    outputFile =fullfile(Dir.results,[subname,'_delay_ft_burst',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.mat']);

    if isfile(matName)
        load(matName);
       [~,chanID] = ismember(chanROI,tfDat{1}.label);

        for cond_i = 1:3
            eeg = tfDat{cond_i};
            parfor ti = 1:height(eeg.trialinfo)
                figure('position',[200 200 800 400]);hold all;axis tight
                imagesc(eeg.time,eeg.freq,squeeze(mean(eeg.powspctrm(ti,chanID,:,:),2)));
                set(gca,'ydir','normal','YTick',[5:5:40],'XTick',[-4.8 -3.6 -2.4 -1.2 0 3],'xlim',SSepochwindow{cond_i})
                xlabel('Time(s)');
                ylabel('Freqeuncy');
                rectangle('Position',[0 0.5 3 40],'LineWidth',1.2,'LineStyle','--')

                hc=colorbar;
                hc.Label.String = 'SD relative to baseline';
                caxis([0 15])
                title([subname, PlotSubj.group{sn},' ',condStr{cond_i}],[chanROI{:}])

                if ~isfolder(fullfile(Dir.figs,subname))
                    mkdir(fullfile(Dir.figs,subname))
                end
                saveas(gca,fullfile(Dir.figs,subname,[subname,condStr{cond_i},'_',[chanROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},num2str(ti),'.png']))
            end
        end
    end
end

