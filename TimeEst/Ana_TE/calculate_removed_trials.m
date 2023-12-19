
clear
load('subs.mat');

List = [];
for sn = 1:height(subs)
    if sn==8 || sn==10
        continue
    end
    subname = subs.name{sn};
    
    set_name = fullfile(Dir.prepro,[subname,'_postICA.set']);
    EEG = pop_loadset('filename',set_name);
    
    %pop_eegplot(EEG,1,1,1);% plot the EEG time series.
    
    List = [List (EEG.trials/600)*100];
    
end
disp(List)
disp(mean(List,"all"))
disp(std(List))