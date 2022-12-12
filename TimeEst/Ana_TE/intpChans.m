function intpChans(sn)

%%
load('subs.mat')

threshV = 8;% for detect bad chans
subname = subs.name{sn};

set_name = fullfile(Dir.prepro,[subname,'_epo.set']);
EEG = pop_loadset('filename',set_name);

[~,tmp1] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',threshV,'norm','on','measure','kurt');
[~,tmp2] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',threshV,'norm','on','measure','prob');


badChanFile = fullfile(Dir.ana,'BadChans',[subname,'_badchans.mat']);
if isfile(badChanFile)
    load(badChanFile)
    disp(badChan)
end

labels = {EEG.chanlocs.labels};
badChan.name = subname;
tmp_chanID = union(tmp1,tmp2);
badChan.autoLabels = labels(tmp_chanID);
fprintf('\rdetected bad channel: %s\r', badChan.autoLabels{:})

%% if there is manually identified bad channels
pop_eegplot(EEG,1,1,1);
tmp = input('badCHans separated by comma: ','s');
badChan.manual = strsplit(tmp,',');
badChan.manual = cellfun(@(x)erase(x," "),badChan.manual,'UniformOutput',false);% remove space if any
if ~isempty(tmp)
    [lia,tmpID] = ismember(badChan.manual,labels);
    bad2interp = unique([tmp_chanID tmpID(lia)']);
else
    bad2interp = tmp_chanID;
end

%%
[~,EOGchanID] = ismember({'TVEOG','BVEOG','LHEOG','RHEOG'},labels);

if ~isempty(bad2interp)
    bad2interp(ismember(bad2interp,EOGchanID))=[];% don't interpolate EOGs
    fprintf('\ninterpolating %s\n',labels{bad2interp});
    EEG = pop_interp(EEG, bad2interp, 'spherical');
    badChan.interped = labels(bad2interp);
end

new_name = [subname,'_pre.set'];
pop_saveset(EEG,'filename',new_name,'filepath',Dir.prepro);
save(badChanFile,'badChan')

