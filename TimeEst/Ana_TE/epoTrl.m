function epoTrl(sn)

%%
load('subs.mat')

subname = subs.name{sn};
set_name = fullfile(Dir.prepro,[subname,'_pre.set']);
EEG = pop_loadset('filename',set_name);

%% time locked to delay onset

markers = {'S 12'}; % cue onset

%     Deci.DT.Markers    = {[28 29 30]};
% Deci.DT.Locks      = [12 17 27];

timelimits = [-1.2 6];

[sEEG,indices] = pop_epoch(EEG,markers,timelimits);
trialN = 600;
rmvd = setdiff(1:trialN,indices);
%% Add new events to a loaded dataset 0.1 second after time-locking event
trialID = 1:trialN;
trialID(rmvd) = [];
if sEEG.trials~=length(trialID)
    error('trial num mismatch')
end

sEEG = eeg_checkset(sEEG, 'eventconsistency');

nevents = length(sEEG.event);
t = 0;
for evt_i = 1 : nevents
    if ischar(sEEG.event(evt_i).type) && strcmpi(sEEG.event(evt_i).type, markers)  && sEEG.event(evt_i).epoch>t

        t = t+1;
        addmarkerStr = ['T' num2str(trialID(t))];
        sEEG.event(end+1) = sEEG.event(evt_i); % Add event to end of event list
        % Specifying the event latency to be 0.1 sec after the referent event (in real data points)
        sEEG.event(end).latency = sEEG.event(evt_i).latency + 0.2*sEEG.srate;
        sEEG.event(end).type = addmarkerStr; % Make the type of the new event
    end
end

[~,I] = sortrows([[sEEG.event.epoch]; [sEEG.event.latency]]',[ 1 2]);
sEEG.event = sEEG.event(I);
sEEG = eeg_checkset(sEEG, 'eventconsistency'); % Check all events for consistency

%%

new_name = [subname,'_epo.set'];
pop_saveset(sEEG,'filename',new_name,'filepath',Dir.prepro);
