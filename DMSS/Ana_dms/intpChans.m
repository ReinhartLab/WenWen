function intpChans(sn)

%%
load('subs.mat')
if subs.rawEEG(sn)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_rawFiltered.set']);

    EEG = pop_loadset('filename',set_name);

    load badChan
    [~,removed_channels] = clean_channels(EEG);

    badChan.name{sn} = subname;
    badChan.labels{sn} = {EEG.chanlocs(removed_channels).labels};
    badChan.labelsID{sn} = find(removed_channels);

    [~,tmp1] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','kurt');
    [~,tmp2] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','prob');
    badChan.indelec{sn} = union(tmp1,tmp2);
    badChan.indelecname{sn} = {EEG.chanlocs(badChan.indelec{sn}).labels};
    save badChan badChan

    bad2interp = union(badChan.labelsID{sn},badChan.indelec{sn}');

    %% if there is manually identified bad channels

    labels = {EEG.chanlocs.labels};
    [~,EOGchanID] = ismember({'TVEOG','BVEOG','LHEOG','RHEOG'},labels);

    if ismember(subname,badChan.name(~cellfun(@isempty,badChan.name)))
        if strcmp(badChan.name{sn},subname)
            Badnames = badChan.chan{sn};
            if ~isempty(Badnames) && ~isempty(Badnames{1})
                [~,badID] =ismember(Badnames,labels);
                bad2interp = union(bad2interp,badID);
            end
        else
            error('mismatched subID and subname')
        end
    end

    if ~isempty(bad2interp)
        bad2interp(ismember(bad2interp,EOGchanID))=[];% don't interpolate EOGs
        EEG = pop_interp(EEG, bad2interp, 'spherical');
        badChan.union{sn} = bad2interp;
        save badChan badChan
    end

    new_name = [subname,'_pre.set'];
    pop_saveset(EEG,'filename',new_name,'filepath',Dir.prepro);

end