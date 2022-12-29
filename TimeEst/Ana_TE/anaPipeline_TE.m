restoredefaultpath % for John's computer

addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))

addpath('D:\intWM-E\toolbox\fieldtrip-20220707')
ft_defaults;
%%
ppts %set project path, read subject folder

%%
clear
load('subs.mat');

parfor sn = 1:height(subs)
    if subs.excluded(sn)==1
        continue
    end
%      loadRaw(sn) 
%      epoTrl(sn)
%      intpChans(sn)
% 
    % rmvArtifact % visual inspection
% 
%     doICA(sn)

    % postICA %
     genCleanData(sn)% will generate new data set
     sanityCheck(sn)% ERP and time freq
end

%%
clear
load('subs.mat');
IsOverWrite = 0;
for IsdePhase = [1 0]
    parfor sn = 1:height(subs)
%         singlePLVsource(sn,IsdePhase,IsOverWrite)

        for IsLap = [0 1]
%             singleTFpower(sn,IsLap,IsdePhase,IsOverWrite)
            singleConn(sn,IsLap,IsdePhase,IsOverWrite)
        end
    end
end

%%

Gplv_ft_sensor

indvPLVpeak% based on singleConn 
indvPLVpeak0.5Hz
Gplv_ft_source %singlePLVsource


