% select the first 5 younger and older participants to estimate sample size

clear;
load('subs.mat');
PilotNames = {'dmss01','dmss02','dmss03','dmss06','dmss04','DMSS030','DMSS031','DMSS032','DMSS036','DMSS037'}';
[~,tmp] = ismember(PilotNames,subs.name);
PilotSubs=subs(tmp,:);
subN = height(PilotSubs);
%%
outA = 3;
RTmax = 5;% RT limit
clear subsBeha

for sub_i = 1:subN
    csvFile = fullfile(Dir.beha,PilotSubs.csvFile{sub_i});

    M = readtable(csvFile);    
    M = M(:,["block_num","ss_num","button_resp_rt","button_resp_corr"]);

    M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});
      
    M(isnan(M.ss_num),:) = [];
    if sum(isnan(M.button_resp_rt)) >0
        M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
    end

    tmp_sum = groupsummary(M,"ss_num","mean",["button_resp_corr"]);
    subsBeha(sub_i,:,1) = table2array(tmp_sum(:,["mean_button_resp_corr"]));%sub*set*(acc)

    tmp_rt = datastats(M.button_resp_rt(M.button_resp_corr==1 & M.button_resp_rt<RTmax));
    downLim = tmp_rt.mean-outA*tmp_rt.std;
    upLim = tmp_rt.mean+outA*tmp_rt.std;

    tmp_idx = M.button_resp_rt<upLim &  M.button_resp_rt>downLim &  M.button_resp_corr==1 & M.button_resp_rt<RTmax;

    tmp_sum = groupsummary(M(tmp_idx,:),"ss_num","mean",["button_resp_rt"]);
    subsBeha(sub_i,:,2) = table2array(tmp_sum(:,["mean_button_resp_rt"]));%sub*set*(rt)
end

PilotSubs.s1_acc = subsBeha(:,1,1);
PilotSubs.s2_acc = subsBeha(:,2,1);
PilotSubs.s3_acc = subsBeha(:,3,1);

PilotSubs.s1_rt = subsBeha(:,1,2);
PilotSubs.s2_rt = subsBeha(:,2,2);
PilotSubs.s3_rt = subsBeha(:,3,2);

beha = PilotSubs;
%%
subN = height(PilotSubs);
setStr = {'Load 1','Load 2','Load 4'};
indStr = {'Accuracy','RT','dprime'};
groupStr = {'Young','Old'};

clear X
X(:,1) = reshape(subsBeha(:,:,1),size(subsBeha,1)*size(subsBeha,2),1);
X(:,2) = repmat(PilotSubs.group,size(subsBeha,2),1);% between
X(:,3) = sort(repmat([1;2;3],subN,1));% within
X(:,4) = repmat([1:subN]',3,1);%subject
[SSQs.acc, DFs.acc, MSQs.acc, Fs.acc, Ps.acc,sst]=mixed_between_within_anova(X);
 
etasq = SSQs.acc{4}/sst;
f = sqrt( etasq / ( 1 - etasq ) )
%%

figure('Name','Beha','Position',[200 200 300 300]);
hold all;axis square;
for g = 1:2
    tmpID = PilotSubs.group==g;
    dat= squeeze(subsBeha(tmpID,:,1));
    mn(g,:) = mean(dat);
    se(g,:) = std(dat)./sqrt(sum(tmpID));    
end

b = bar(mn','BarWidth',1);
errorbar(vertcat(b.XEndPoints)',mn',se','k.','LineWidth',1)
set(gca,'XTickLabel',setStr,'XLim',[0.25 3.75],'XTick',1:3)

legend(groupStr,'Location','southoutside')
ytickformat('%,.1f')