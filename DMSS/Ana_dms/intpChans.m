function intpChans(sn)

%%
load('subs.mat')
if subs.rawEEG(sn)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_rawFiltered.set']);

    EEG = pop_loadset('filename',set_name);

    matName = fullfile(Dir.ana,'preBadChans',[subname,'_preInterp.mat']);
    if isfile(matName)
        load(matName)
    end
    [~,removed_channels] = clean_channels(EEG);

    badChan.labels = {EEG.chanlocs(removed_channels).labels};
    badChan.labelsID = find(removed_channels);

    [~,tmp1] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','kurt');
    [~,tmp2] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','prob');
    badChan.indelec = union(tmp1,tmp2);
    badChan.indelecname = {EEG.chanlocs(badChan.indelec).labels};
    save(matName,'badChan')

    bad2interp = union(badChan.labelsID,badChan.indelec');

    %% if there is manually identified bad channels

    labels = {EEG.chanlocs.labels};
    [~,EOGchanID] = ismember({'TVEOG','BVEOG','LHEOG','RHEOG'},labels);

    if strcmp(badChan.name,subname)
        if isfield(badChan,'chan')
            Badnames = badChan.chan;
            if ~isempty(Badnames) && ~isempty(Badnames{1})
                [~,badID] =ismember(Badnames,labels);
                bad2interp = union(bad2interp,badID);
            end
        end
    else
        error('mismatched subID and subname')
    end

    if ~isempty(bad2interp)
        bad2interp(ismember(bad2interp,EOGchanID))=[];% don't interpolate EOGs
        EEG = pop_interp(EEG, bad2interp, 'spherical');
        badChan.union = bad2interp;
        save(matName,'badChan')
    end

    new_name = [subname,'_pre.set'];
    pop_saveset(EEG,'filename',new_name,'filepath',Dir.prepro);

end