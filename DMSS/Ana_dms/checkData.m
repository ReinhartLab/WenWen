function checkData(sn)
threshV = 100;
locthresh = 8;
globthresh = 8;
zthresh = 8;

load('subs.mat');

if subs.rawEEG(sn)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_postICA.set']);
    EEG = pop_loadset('filename',set_name);

    pop_eegplot(EEG,1,1,1);
    %% interpolate bad channels after ICA
    [~,tmp1] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','kurt');
    [~,tmp2] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','prob');
    tmpAutoBadID = union(tmp1,tmp2);
    labels = {EEG.chanlocs.labels};
    badChan2interp.badChansAuto = labels(tmpAutoBadID);

    postIterpFile = fullfile(Dir.ana,'postICAinterpChans',[subname,'_postInterp.mat']);

    if isfile(postIterpFile)
        load(postIterpFile)
        if isfield(badChan2interp,'manual')
            disp(badChan2interp.manual)
        end
    end
    tmp = input('badCHans separated by comma: ','s');
    badChan2interp.manual = strsplit(tmp,',');

    if ~isempty(tmp)
        [~,tmpID] = ismember(badChan2interp.manual',labels);
        tmpbad =  unique([tmpAutoBadID tmpID]);
    else
        tmpbad =  tmpAutoBadID;
    end
    
    [~,EOGchanID] = ismember({'TVEOG','BVEOG','LHEOG','RHEOG'},labels);
    tmpbad(ismember(tmpbad,EOGchanID))=[];% don't interpolate EOGs
    badChan2interp.interp = tmpbad;

    disp(labels(tmpbad))

    EEG = pop_interp(EEG, tmpbad, 'spherical');
    save(postIterpFile,'badChan2interp')
    %% rejection
    disp('Removing EOGs before trials rejection')
    sEEG = pop_select(EEG,'time',[-5 5],'nochannel',{'TVEOG','LHEOG','RHEOG','BVEOG'});% remove EOG channels

    [~, rejT] = pop_eegthresh( sEEG, 1, 1:sEEG.nbchan, -threshV, threshV, sEEG.xmin, sEEG.xmax,0,0);
    [~, ~, ~, rejK] = pop_rejkurt_ww( sEEG, 1, 1:sEEG.nbchan, locthresh,globthresh, 'reject',0,'plotflag',0);
    [~, ~, ~, rejP] = pop_jointprob_ww( sEEG, 1, 1:sEEG.nbchan, locthresh, globthresh,1, 'reject',0,'plotflag',0);

    tmp = abs(zscore(sEEG.data,0,3))>zthresh;% along trial dimension
    rejZ = squeeze(sum(sum(tmp),2))>0;

    tmp = abs(zscore(sEEG.data,0,2))>zthresh;% along time dimension
    rejZ2 = squeeze(sum(sum(tmp),2))>0;

    %% summary
    rej_union = union(rejT,[find(rejK) find(rejP) find(rejZ)' find(rejZ2)']);
    fprintf(['\n------trials summary-----\n' ...
        'Rejecting %d trials of %d\n' ...
        'Remaining %d of 720 (%.3f)\n'], ...
        length(rej_union),EEG.trials,EEG.trials-length(rej_union),(EEG.trials-length(rej_union))/720);

    %% remove

    EEG = pop_rejepoch( EEG, rej_union,0);
    pop_eegplot(EEG,1,1,1);

    %% save

    set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
    pop_saveset(EEG,'filename',set_name);

end
end

