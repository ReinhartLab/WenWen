function sanityCheck(sn)
load('subs.mat');
subname = subs.name{sn};
if subs.excluded(sn)==1
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
EEG = pop_loadset('filename',set_name);

EEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

%%  re-aligned to feedback onset:29-early; 30-late; 28-correct

marker ={'S 28','S 29','S 30'};
EEG = pop_epoch(EEG,marker,[-1 2]);

%% plot
% topoT  = input('enter time samples(in ms) with brackets[]:');
topoT = [300 500];
figure('position',[200 200 1000 600]);
subplot(1,2,1)
pop_timtopo(EEG, 1000*[EEG.xmin EEG.xmax], topoT);

subplot(1,2,2)
pop_spectopo(EEG,1,[0.1 0.6]*1000,'ERP','freq',[2:2:10],'freqrange',[1 40],'title',['sn#',num2str(sn)])
saveas(gca,fullfile(Dir.figs,'SanityCheck',sprintf('sn%02d_%s.png',sn,subname)))