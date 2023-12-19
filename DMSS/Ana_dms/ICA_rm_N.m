clear
load('subs.mat');

trialN = nan(height(subs),1);
origN = nan(height(subs),1);
for sn = 1:height(subs)

    subname = subs.name{sn};
    rejICA = fullfile(Dir.ana,'rmvICA',[subname,'_rmvID.mat']);

    if isfile(rejICA) && subs.exclude(sn)==0
        load(rejICA)
        subs.ICN(sn) = length(rmvd);
    end
end

grpstats(subs(subs.exclude==0 & subs.rawEEG ==1,:),'groupStr',["mean","std","range","max","min"],'DataVars',"ICN")

 [~,p,~,s] = ttest2(subs.ICN(subs.exclude==0 & subs.rawEEG ==1 & subs.group ==1),subs.ICN(subs.exclude==0 & subs.rawEEG ==1 & subs.group ==2))