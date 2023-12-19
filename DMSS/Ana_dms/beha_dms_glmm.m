clear;
load subs.mat
subs(subs.rawEEG==0 | subs.exclude ==1,:) = [];

subN = height(subs);
IsExcludeOutlierRT = 0;
%%
outA = 3;
RTmax = 5;% RT limit
myTbl = table;

for sub_i = 1:subN
    csvFile = fullfile(Dir.beha,subs.csvFile{sub_i});

    M = readtable(csvFile);    
    M = M(:,["block_num","ss_num","button_resp_rt","button_resp_corr"]);

    M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});
      
    M(isnan(M.ss_num),:) = [];
    if sum(isnan(M.button_resp_rt)) >0
        M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
    end

    if IsExcludeOutlierRT ==1         %---exclude trials
    tmp_rt = datastats(M.button_resp_rt(M.button_resp_corr==1 & M.button_resp_rt<RTmax));
    downLim = tmp_rt.mean-outA*tmp_rt.std;
    upLim = tmp_rt.mean+outA*tmp_rt.std;
    tmp_idx = M.button_resp_rt<upLim &  M.button_resp_rt>downLim &  M.button_resp_corr==1 & M.button_resp_rt<RTmax;
    M(tmp_idx,:)= [];
    end

    M.subj = repmat(subs.name(sub_i),height(M),1);
    M.group = repmat(subs.groupStr(sub_i),height(M),1);
    myTbl = [myTbl;M];   
end

%%
myTbl.ss_num = categorical(myTbl.ss_num);
glme = fitglme(myTbl,'button_resp_corr ~ ss_num*group +(1|subj)','Distribution','Binomial');
glme.anova

