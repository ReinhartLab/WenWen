function intpChans(sn)

%%
load('subs.mat')

subname = subs.name{sn};
set_name = fullfile(Dir.prepro,[subname,'_rawFiltered.set']);
EEG = pop_loadset('filename',set_name);

[~,removed_channels] = clean_channels(EEG);
[~,tmp1] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','kurt');
[~,tmp2] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',8,'norm','on','measure','prob');

labels = {EEG.chanlocs.labels};

badChan.name = subname;
badChan.chanID = union(tmp1,[tmp2 find(removed_channels)]);
badChan.labels = labels(badChan.chanID);
fprintf('\rbad channels detected: %s\r', badChan.labels{:})
%% if there is manually identified bad channels

tmp = input('badCHans separated by comma: ','s');
badChan.manual = strsplit(tmp,',');
if ~isempty(tmp)
    [~,tmpID] = ismember(badChan.manual,labels);
    bad2interp = [badChan.chanID,tmpID];
else
    bad2interp = badChan.chanID;
end

badChanFile = fullfile(Dir.ana,'BadChans',[subname,'badchans.mat']);
if isfile(badChanFile)
    load(badChanFile)
    disp(badChan)
end
%%
[~,EOGchanID] = ismember({'TVEOG','BVEOG','LHEOG','RHEOG'},labels);

if ~isempty(bad2interp)
    bad2interp(ismember(bad2interp,EOGchanID))=[];% don't interpolate EOGs
    EEG = pop_interp(EEG, bad2interp, 'spherical');
    badChan.bad2interp = bad2interp;
end

new_name = [subname,'_pre.set'];
pop_saveset(EEG,'filename',new_name,'filepath',Dir.prepro);
save(badChanFile,'badChan')

