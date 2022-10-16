% visual inspection bad channels, be conservative
% this will be combined with auto bad channel rejection
clear
close all
load('subs.mat');
load badChan
%%
for sn = 29%1:height(subs)
    if subs.rawEEG(sn)==1
        subname = subs.name{sn};

        %%
        set_name = fullfile(Dir.prepro,[subname,'_rawFiltered.set']);
        EEG = pop_loadset('filename',set_name);

        sEEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

        pop_eegplot(sEEG,1,1,1);
        %%
        badChan.name{sn}=subname;
        disp(badChan.chan{sn})
        tmp = input('badCHans separated by comma: ','s');
        if ~isempty(tmp)
            badChan.chan{sn} = strsplit(tmp,',');
        end
        save badChan badChan

    end
end


