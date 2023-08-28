function sanityCheck(sn)
load('subs.mat');
subname = subs.name{sn};
if subs.exclude(sn)==1 || subs.rawEEG(sn)==0
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
EEG = pop_loadset('filename',set_name);

EEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

%%  re-aligned to delay onset

if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
    stim1 = {'S 11','S 21', 'S 31'};
else
    stim1 = {'11','21','31'};
end
%% load Beha
csvFile = fullfile(Dir.beha,subs.csvFile{sn});
M = readtable(csvFile);
M = M(:,["block_num","ss_num","button_resp_rt","button_resp_corr"]);
%                                M = M(:,["block_num","ss_num","type","button_resp_rt","button_resp_corr"]);
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
figure('position',[200 200 1800 900]);
ssN = [1 2 4];

toi = {[-1.4 3.5],[-2.6 3.5],[-5 3.5]};% time-locked2delay
bsWindow = {[-1.4 -1.2],[-2.6 -2.4],[-5 -4.8]};
for ss = 1:3
    sEEG = pop_select(EEG,'trial',find(M.ss_num == ssN(ss)),'time',toi{ss});

    sEEG = pop_rmbase(sEEG,bsWindow{ss}*1000);% pre trial baseline

    % topoT  = input('enter time samples(in ms) with brackets[]:');
    topoT = [500 1500 2500];

    subplot(2,3,ss)
    pop_timtopo(sEEG, 1000*[sEEG.xmin sEEG.xmax], topoT);

    subplot(2,3,ss+3)
    pop_spectopo(sEEG,1,[0.1 0.6]*1000,'ERP','freq',[5:5:25],'freqrange',[1 40],'title',['sn',num2str(sn),subname])

end
saveas(gca,fullfile(Dir.figs,'SanityCheck',sprintf('sn%02d_%s.png',sn,subname)))