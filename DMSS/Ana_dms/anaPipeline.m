%% set up toolbox
restoredefaultpath % for John's computer

addpath(genpath('D:\intWM-E\toolbox\eeglab2021.1'))
addpath('D:\intWM-E\toolbox\fieldtrip-20220707')

ft_defaults;
%%
ppts %set project path, read subject folder

beha_dms % behavior

%% preprocessing
clear
load('subs.mat');

for sn = 1:height(subs)
    %     loadRaw(sn) %
    % badChanInspection % manually reject bad channels

    %     intpChans(sn)% automatic
    %     epoTrl(sn)

    % rmvArtifact % visual inspection

    %         doICA(sn)

    % postICA
    %     checkData(sn)% will generate new data set
    %     sanityCheck(sn) %gen figures
end

% checkEEGMarkerBehaMatching
% cleanTrialsCalculator % will add a column in subs

%% run analysis at individual level

clear
load('subs.mat');
IsOverwrite = 0; % overwrite existing output?
for IsCorretTrials = [1 0] % using correct or all trials for analysis
    for IsdePhase = [1 0] % whether to subtract ERP from each trial for time freq decomposition
        M = 5;
        parfor (sn = 1:height(subs),M)
            for IsBL2preDelay = [0 1] % baselined to preTrial interval or preDelay/preResponse interval

                singleVariance_delay_ft_LCMV(sn,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                singleVariance_response_ft_LCMV(sn,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                for IsLap = [1 0] % laplacian filter
                    single_variablity_delay_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                    singleTF_delay_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                    singleTF_response_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                    single_variablity_response_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                    single_variablity_response_ft_oddeven(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                    iterN = 200;
                    single_variablity_response_ft_correctVSwrong(sn,IsLap,IsdePhase,0,IsBL2preDelay,IsOverwrite,iterN)
                    single_power_response_ft_correctVSwrong(sn,IsLap,IsdePhase,0,IsBL2preDelay,IsOverwrite,iterN)

                    %                     singleTF_delay_burst_ft(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                    %                     singleTF_delay_ITPC_ft_wavelet(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)

                    single_power_delay_ft_trialCorrelation_Nsampled(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)% to do
                    for Nsamps = [6:2:12]% size of trial cluster
                        single_variablity_delay_ft_trialCorrelation_Nsampled(sn,Nsamps,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                    end

                    single_power_response_ft_Nplus1Correlation_Nsampled(sn,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)% to do
                    for Nsamps = [12:4:30]% size of trial cluster, odd even alternated
                        single_variablity_response_ft_Nplus1Correlation_Nsampled(sn,Nsamps,IsLap,IsdePhase,IsCorretTrials,IsBL2preDelay,IsOverwrite)
                    end
                end
            end
        end
    end
end


%% decoding

clear
load('subs.mat');

nIter = 50;
nfolds = 10;
IsOverwrite = 0;
for IsOcci = [0 1]
    for IsBL2preDelay = [0 1]
        for IsCorretTrials = [1 0]
            parfor sn = 1:height(subs)

                singleSSdelay_mal(sn,nIter,nfolds,IsOcci,IsBL2preDelay,IsCorretTrials,IsOverwrite)
                singleSSdelaySearchLight_mal(sn,nIter,nfolds,IsBL2preDelay,IsCorretTrials,IsOverwrite)

                for IsdePhase = [1 0]
                    singleSSdelayFreq_mal(sn,nIter,nfolds,IsOcci,IsdePhase,IsBL2preDelay,IsCorretTrials,IsOverwrite)
                end
            end
        end
    end
end
%% plotting

GndStimERP

GndSSdelay
GndSSdelaySearchLight
GndSSdelayFreq

GndNeuralVar
GndNeuralVar_encoding
GndNeuralVar_response
GndNeuralVar_responseOddEven
GndNeuralVar_responseCorVSwrong % iterN = 200
GndNeuralVar_response_Nsampled
GndNeuralVar_delay_Nsampled_correctTrialsRT % using correct trials, behavioral index = RT since ACC would be one
GndNeuralVar_delay_Nsampled_allTrials % include error trials

GndTF_ft
GndTF_ft_preDelayBL % IsBL2preDelay =1, different Ylim
GndTF_encoding
GndTF_ft_response
GndTF_responseCorVSwrong 
GndTF_response_Nsampled
GndTF_delay_Nsampled

GndTF_ITPC_ft
GndPlotSingleTrialBurst

GndPlotSource
GndPlotSource_Resp

