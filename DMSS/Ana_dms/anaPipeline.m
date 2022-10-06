restoredefaultpath % for John's computer

addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))

addpath('D:\intWM-E\toolbox\fieldtrip-20211209')

%%
% ppts %set project path, read subject folder, input bad channels for EEG interpolate

% beha_dms

badChanInspection %
%%
clear
load('subs.mat');

parfor sn = 1:height(subs)
%         loadRaw(sn)
% 
%     intpChans(sn)
%      epoTrl(sn)

    % rmvArtifact % visual inspection

    doICA(sn)

    % postICA %
     checkData(sn)% will generate new data set

end
% reEpoch2delay
%%
% clear
% load('subs.mat');
% 
% for IsLap = [1 0]
%     for IsdePhase = [1 0]
%         for IsCorretTrials = [1 0]
%             for IsBL2preDelay = [0 1]
%                 parfor sn = 1:height(subs)
%                     singleTF_delay_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay)
% 
%                 end
%             end
%         end
%     end
% end


%%

% clear
% load('subs.mat');
% 
% nIter = 50;
% nfolds = 10;
%  
% for IsBL2preDelay = [1 0]
%     for IsOcci = [1 0]
%         for IsCorretTrials = [1 0]
%             parfor sn = 1:height(subs)
% %                 singleSSdelay_mal(sn,nIter,nfolds,IsOcci,IsBL2preDelay,IsCorretTrials)
%                 singleSSdelaySearchLight_mal(sn,nIter,nfolds,IsBL2preDelay,IsCorretTrials)
% % 
% %                 for IsdePhase = [1 0]
% %                     singleSSdelayFreq_mal(sn,nIter,nfolds,IsOcci,IsdePhase,IsBL2preDelay,IsCorretTrials)
% %                 end
%             end
%         end
%     end
% end
%%





