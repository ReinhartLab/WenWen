clear;rng('shuffle');
load('subs.mat');

mySubs = 1:subs.subN;
mySubs([5]) = [];
subN = length(mySubs);

myFigBasic
addpath(genpath('D:\Toolbox\fieldtrip-20211209'))
ft_defaults;

addpath('D:\OneDrive - Boston University\PostDoc\Research\WMinterruption\ana_intWM\f_PlotEEG_BrainNetwork')

IsLap = 1;
IsdePhase = 1;
IsCorretTrials = 1;
txtCell = {'','','';'_lap','_dephase','_corrTrials'};
permFile = ['Perm_dotsPLV',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat'];

%%
clear tfAll gndTF
for sub_i = 1:subN
    sn = mySubs(sub_i);
    subname = subs.name{sn};

    load(fullfile(Dir.results,[subname,'_coh',txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},'.mat']))

    for cond_i = 1:3
        dat = cat(4,tfDat{cond_i}.plv.data{:});%time*chan*chan*freq;
        times = tfDat{1}.plv.times;
        BLtime = [-200 0];
        BLidx = dsearchn(times',BLtime');
        BLidx = BLidx(1):BLidx(2);
        dat = bsxfun(@minus,dat,mean(dat(BLidx,:,:,:)));

        tfAll.subs{cond_i}(sub_i,:,:,:,:) = zscore(dat,1);
    end
end

gndTF = cellfun(@(x)squeeze(mean(x)),tfAll.subs,'UniformOutput',false);

condStr = {'Fixation','Ignore','Attend'};
condDiffStr = {'Ignore-Fixation','Attend-Fixation','Attend-Ignore','Ignore-Attend'};

%%
myFigBasic
clear chansID
chanLabels = tfDat{1}.itc.label;
% chans{1} = {'FC4','F4'};
chans{1} = {'FC3','F3'};

chans{2} = {'POz','POz'};
% chans{1} = {'Fz','Fz','FCz','FCz'};
% chans{2} = {'POz','Pz','POz','Pz'};
[log1,chansID(:,1)] = ismember(chans{1},chanLabels);
[log2,chansID(:,2)] = ismember(chans{2},chanLabels);
chansID = [chansID(log1,1) chansID(log2,2)];
chansID = sort(chansID,2);% The smaller channel number is to be used first.

freqs = tfDat{1}.plv.freq(2,:);
figure;
for cond = 1:3
    subplot(3,1,cond);
    clear dat
    for c = 1:size(chansID,1)
        dat(c,:,:) = squeeze(gndTF{cond}(:,chansID(c,1),chansID(c,2),:))';%freq* time
    end

    dat = squeeze(mean(dat));
    imagesc(times,freqs,dat);
    rectangle('Position',[0 0.5 2000 40],'LineWidth',1.2,'LineStyle','--')
    set(gca,'YDir','normal','YTick',freqs(1:3:end),'TickDir','in','TickLength',[0.01 0.01])
    xlabel('\bfTime(ms)')
    ylabel('\bfFrequency(Hz)')

%     title([condStr{cond} ':FC4&F4-POz' ])
        title([condStr{cond} ':FC3&F3-POz' ])

    cb = colorbar;
    cb.Label.String = 'Normalized PLV';
    cb.Limits = [-0.8 0.8];
end

%%
toi = [2 2.5]*1000;
toi = times<=toi(2) & times>=toi(1);
% myFreqID = [15 25];%Hz

myFreqID = [9 13];%Hz
% myFreqID = [3 7];%Hz

myFreqID = mean(tfDat{1}.plv.freq)<=myFreqID(2) & mean(tfDat{1}.plv.freq)>=myFreqID(1);
figure
for cond_i = 1:3

    aij = squeeze(mean(mean(gndTF{cond_i}(toi,:,:,myFreqID),1),4));

    nch = 58; %take care of this variable (nch must be according to matrix size you want to plot it)
    p = 0.01;   %proportion of weigthed links to keep for.
    aij = threshold_proportional(aij, p); %thresholding networks due to proportion p
    ijw = adj2edgeL(triu(aij));             %passing from matrix form to edge list form
    n_features = sum(aij, 2);
    subplot(1,3,cond_i)
    f_PlotEEG_BrainNetwork(nch, ijw, 'w_wn2wx', n_features, 'n_nn2nx', 'ncb');
    title([condStr{cond_i} '-alpha'])
%             title([condStr{cond_i} '-theta'])

end
%%

for cond_i = 1:3
    aijDat{cond_i} = squeeze(mean(mean(gndTF{cond_i}(toi,:,:,myFreqID),1),4));
end
aijDiff{1} = aijDat{2}-aijDat{1};
aijDiff{2} = aijDat{3}-aijDat{1};
aijDiff{3} = aijDat{3}-aijDat{2};
aijDiff{4} = aijDat{2}-aijDat{3};

figure
cond2plot = [1 2 3];
% cond2plot = [1 2 4];

% cond2plot = [1 2 3 4];

for cond_i = 1:length(cond2plot)

    tmpC = cond2plot(cond_i);
    aij=aijDiff{tmpC};

    nch = 58; %take care of this variable (nch must be according to matrix size you want to plot it)
    p = 0.01;   %proportion of weigthed links to keep for.
    aij = threshold_proportional(aij, p); %thresholding networks due to proportion p
    ijw = adj2edgeL(triu(aij));             %passing from matrix form to edge list form
    n_features = sum(aij, 2);

    subplot(1,length(cond2plot),cond_i)

    f_PlotEEG_BrainNetwork(nch, ijw, 'w_wn2wx', n_features, 'n_nn2nx', 'ncb');
%         title([condDiffStr{cond_i} '-theta'])
    title([condDiffStr{tmpC} '-alpha'])
end


