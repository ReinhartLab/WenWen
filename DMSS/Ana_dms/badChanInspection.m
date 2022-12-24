% visual inspection bad channels, be conservative
% this will be combined with auto bad channel rejection
clear
close all
load('subs.mat');
%%
for sn = 21%1:height(subs)
    if subs.rawEEG(sn)==1
        subname = subs.name{sn};

        %%
        set_name = fullfile(Dir.prepro,[subname,'_rawFiltered.set']);
        EEG = pop_loadset('filename',set_name);

        sEEG = pop_select(EEG,'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

        pop_eegplot(sEEG,1,1,1);
        %%
        matName = fullfile(Dir.ana,'preBadChans',[subname,'_preInterp.mat']);
        if isfile(matName)
            load(matName)
            if istable(badChan)
                badChan = table2struct(badChan);
            end
            if isfield(badChan,'chan')
                disp(badChan.chan)
            end
        end
        badChan.name=subname;
        tmp = input('badCHans separated by comma: ','s');
        if ~isempty(tmp)
            badChan.chan = strsplit(tmp,',');
        end
        save(matName,'badChan')
    end
end


