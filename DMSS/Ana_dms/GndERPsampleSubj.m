clear;
load('subs.mat');
load beha.mat
IsCorrect = 1;
ssN = [1 2 4];

txtCell = {'','';'_corrTrials','_bl2preDelay'};
load(['ERP',txtCell{IsCorrect+1,1},'.mat'])

%%
subs(subs.rawEEG==0,:) = [];
subsBeha(subs.rawEEG==0,:,:) = [];
subsAll(cellfun(@isempty,subsAll(:,1)),:) = [];
%%
iter = 0;
nIter = 500;
forT = nan(nIter,1);
pickSubsAll = nan(nIter,sum(subs.group==2));
while iter<nIter
    pickYoung = randperm(sum(subs.group==1),sum(subs.group==2));
    youngBeha = subsBeha(pickYoung+sum(subs.group==2),:,:);
    oldBeha = subsBeha(subs.group==2,:,:);

    [h,p,~,stat] = ttest2(youngBeha,oldBeha);
    if sum(h)==0
        iter = iter+1;
        pickSubsAll(iter,:) = pickYoung;
        forT(iter) = sum(abs(stat.tstat),"all");
    end
end
[~,tmp] = min(forT);
pickYoung = pickSubsAll(tmp,:);
excludeYoung = 1:sum(subs.group==1);
excludeYoung(pickYoung) = [];
excludeYoung = excludeYoung+sum(subs.group==2);
%%

subs(excludeYoung,:) = [];
subsBeha(excludeYoung,:,:) = [];
subsAll(excludeYoung,:) = [];

%%
