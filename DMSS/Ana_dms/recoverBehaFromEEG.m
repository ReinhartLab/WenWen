clear;
load subs.mat
subname = 'DMSS032';
set_name = fullfile(Dir.prepro,[subname,'_rawFiltered.set']);
EEG = pop_loadset('filename',set_name);


event = struct2table(EEG.event);
% tmp = cellfun(@(x)contains(x,'Stimulus'),[event.code]);
tmp = cellfun(@(x)contains(x,'Stim'),[event.code]);
event = event(tmp,:);

%% loop
fbCode = {'S 72','S 73','S 74'};%load1 2 4
probeCode = {'S 57'};% probe on
respCode = {'S 60'};% key is done
buttonCode = {'S 61','S 62'};% left right button
corrCode = {'S 64','S 65'};% correct incorrect
[tmp,a] = cellfun(@(x)ismember(x,fbCode),[event.type]);
ss = a(tmp);

[tmp,a] = cellfun(@(x)ismember(x,corrCode),[event.type]);
ACC = 2-a(tmp);

[tmp,a] = cellfun(@(x)ismember(x,buttonCode),[event.type]);
button = a(tmp)+3;

tmpProbe = cellfun(@(x)ismember(x,probeCode),[event.type]);
tmpKey = cellfun(@(x)ismember(x,respCode),[event.type]);

RT = event.latency(tmpKey)-event.latency(tmpProbe);

beha = table(RT,ss,button,ACC);