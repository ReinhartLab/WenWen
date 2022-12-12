clear;
close all

load('subs.mat');

sn = 3;
checkData(sn)

%%
subname = subs.name{sn};
set_name = fullfile(Dir.prepro,[subname,'_ICA.set']);
EEG = pop_loadset('filename',set_name);

% EEG = pop_iclabel(EEG,'default');
% pop_selectcomps(EEG, [1:min([45 size(EEG.icaweights,1)])]);
% pop_eegplot( EEG, 0, 1, 1);

rejICA = fullfile(Dir.ana,'rmvICA',[subname,'_rmvID.mat']);
if isfile(rejICA)
    load(rejICA)
    disp(rmvd)
end

pause
rmvd  = input('enter to-be-rmvd cmp ID with brackets[]:');

EEG = pop_subcomp(EEG,rmvd,0); % change to 1 to confirm rejection

set_name = fullfile(Dir.prepro,[subname,'_postICA.set']);
pop_saveset(EEG,'filename',set_name);
save(rejICA,'rmvd')

% checkData(sn)
