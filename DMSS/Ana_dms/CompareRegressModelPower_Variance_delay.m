
clear;rng('shuffle');
load('subs.mat');

IsCorretTrials = 1;
IsBL2preDelay = 0;% baselined to pre response
IsLap = 1;
IsdePhase=1;

txtCell = {'','','','';'_lap','_dephase','_corrTrials','_bl2preDelay'};

groupStr = {'Young','Old','Young-Old'};
condStr = {'ss1','ss2','ss4'};

frontalROI = {'Fz','F1','F2','FCz','AFz'};% frontal cluster
% frontalROI = {'CPz','CP1','CP2','Pz','Cz'};%centroparietal cluster

freq.betaFreq = [15 25];% Hz
timeROI.Post = [0 3];% in s

varTable=readtable(fullfile(Dir.ana,'TableOutput',['DelayBetaVariance_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']));
powTable=readtable(fullfile(Dir.ana,'TableOutput',['DelayBetaPower_',num2str(freq.betaFreq(1)),'~',num2str(freq.betaFreq(2)),'Hz',num2str(timeROI.Post(1)),'~',num2str(timeROI.Post(2)),'s',[frontalROI{:}],txtCell{IsLap+1,1},txtCell{IsdePhase+1,2},txtCell{IsCorretTrials+1,3},txtCell{IsBL2preDelay+1,4},'.csv']));


load beha.mat
wanted = ismember(beha.name,powTable.name);
beha = beha(wanted,:);
behaStr = {'RT','acc'};

addpath(genpath('D:\Toolbox\crameri_v1.08'))
myColors = crameri('bamako',2);% https://www.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps
%%

for idx = 2%1:2
    if idx == 2
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_acc','s2_acc','s3_acc'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_acc','s2_acc','s3_acc'}));
    else
        dat.beha{1} = table2array(beha(beha.group==1,{'s1_rt','s2_rt','s3_rt'}));
        dat.beha{2} = table2array(beha(beha.group==2,{'s1_rt','s2_rt','s3_rt'}));
    end

    for gi = 1%1:2

        clear X Y
        Y = reshape(dat.beha{gi},size(dat.beha{gi},1)*size(dat.beha{gi},2),1);
        X(:,1) = sort(repmat([1 2 4]',size(dat.beha{gi},1),1));
        X(:,2) = [powTable.load1(powTable.group==gi);powTable.load2(powTable.group==gi);powTable.load3(powTable.group==gi)];
        mdl_124_p = fitglm(X,Y);
        
        X(:,2) = [varTable.load1(varTable.group==gi);varTable.load2(varTable.group==gi);varTable.load3(varTable.group==gi)];
        mdl_124_v =  fitglm(X,Y);
        
        mdl_124_p.Coefficients
        mdl_124_p.Rsquared

        mdl_124_v.Coefficients
        mdl_124_v.Rsquared

        X(:,1) = sort(repmat([1 2 4]',size(dat.beha{gi},1),1));
        X(:,2) = [powTable.load1(powTable.group==gi);powTable.load2(powTable.group==gi);powTable.load3(powTable.group==gi)];
        X(:,3) = [varTable.load1(varTable.group==gi);varTable.load2(varTable.group==gi);varTable.load3(varTable.group==gi)];
        mdl_124_s =  stepwiseglm(X,Y,'constant','Distribution','normal','Criterion','adjrsquared','Lower','y ~ x1','upper','linear');

        clear X Y
        Y = reshape(dat.beha{gi}(:,[2 3]),size(dat.beha{gi},1)*2,1);
        X(:,1) = sort(repmat([2 4]',size(dat.beha{gi},1),1));
        X(:,2) = [powTable.load2(powTable.group==gi);powTable.load3(powTable.group==gi)];
        mdl_24_p =  fitglm(X,Y);
     
        X(:,2) = [varTable.load2(varTable.group==gi);varTable.load3(varTable.group==gi)];
        mdl_24_v =  fitglm(X,Y);

        mdl_24_p.Coefficients
        mdl_24_p.Rsquared

        mdl_24_v.Coefficients
        mdl_24_v.Rsquared
    end
end
