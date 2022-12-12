restoredefaultpath % for John's computer

addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))

addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
%%
ppts %set project path, read subject folder

%%
clear
load('subs.mat');

for sn = [13 14]%1:height(subs)
    if subs.excluded(sn)==1
        continue
    end
     loadRaw(sn) %done
     epoTrl(sn)
     intpChans(sn)

    % rmvArtifact % visual inspection

    doICA(sn)

    % postICA %
%      checkData(sn)% will generate new data set

end

%%
clear
load('subs.mat');
IsOverWrite = 0;
for IsLap = [1 0]
    for IsdePhase = [1 0]
        parfor sn = 1:height(subs)
            singleConn(sn,IsLap,IsdePhase,IsOverWrite)
        end
    end
end


%%


%%




