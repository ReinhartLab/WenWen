restoredefaultpath % for John's computer

addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))

addpath('D:\intWM-E\toolbox\fieldtrip-20211209')
ft_defaults;
%%
ppts %set project path, read subject folder

beha_dms

%%
clear
load('subs.mat');

for sn = 1:height(subs)
%     loadRaw(sn)
% badChanInspection % manually reject bad channels

%     intpChans(sn)% automatic
%     epoTrl(sn)

    % rmvArtifact % visual inspection

%         doICA(sn)

    % postICA %
    checkData(sn)% will generate new data set
end

% checkEEGMarkerBehaMatching
% cleanTrialsCalculator % will add a column in subs
%%
clear
load('subs.mat');
IsOverwrite = 1;
for IsLap = [1 0]
    for IsdePhase = [1 0]
        for IsCorretTrials = [1 0]
            parfor sn = 1:height(subs)
                for IsBL2preDelay = [0 1]

                                        singleTF_delay_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                                        singleTF_response_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

%                                         singleTF_delay_ft_DICS(sn,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)


                                        single_variablity_delay_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                    single_variablity_response_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                    singleTF_delay_burst_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                    singleTF_delay_ITPC_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                end
            end
        end
    end
end


%%

clear
load('subs.mat');

nIter = 50;
nfolds = 10;
IsOverwrite = 1;
for IsOcci = [0 1]
    for IsBL2preDelay = [0 1]
        for IsCorretTrials = [1 0]
            parfor sn = 1:height(subs)

                %                 singleSSdelay_mal(sn,nIter,nfolds,IsOcci,IsBL2preDelay,IsCorretTrials,IsOverwrite)
                singleSSdelaySearchLight_mal(sn,nIter,nfolds,IsBL2preDelay,IsCorretTrials,IsOverwrite)

                for IsdePhase = [1 0]
                    singleSSdelayFreq_mal(sn,nIter,nfolds,IsOcci,IsdePhase,IsBL2preDelay,IsCorretTrials,IsOverwrite)
                end
            end
        end
    end
end
%%

GndStimERP

GndSSdelay
GndSSdelaySearchLight
GndSSdelayFreq

GndNeuralVar
GndNeuralVar_response

GndTF_ft
GndTF_ft_response

GndTF_ITPC_ft
GndPlotSingleTrialBurst


