restoredefaultpath % for John's computer

addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))

addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
%%
ppts %set project path, read subject folder

%%
clear
load('subs.mat');

parfor sn = 1:height(subs)
     loadRaw(sn)
%      intpChans(sn)
%      epoTrl(sn)

    % rmvArtifact % visual inspection

%     doICA(sn)

    % postICA %
%      checkData(sn)% will generate new data set
end

%%
clear
load('subs.mat');

for IsLap = [1 0]
    for IsdePhase = [1 0]
        for IsCorretTrials = [1 0]
            for IsBL2preDelay = [0 1]
                parfor sn = 1:height(subs)
                    singleTF_delay_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay)
% single_connectivity_delay_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay)
                end
            end
        end
    end
end


%%


%%

GndStimERP

GndSSdelay
GndSSdelaySearchLight
GndSSdelayFreq



