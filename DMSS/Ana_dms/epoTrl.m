function epoTrl(sn)

%%
load('subs.mat')
if subs.rawEEG(sn)

    subname = subs.name{sn};
    set_name = fullfile(Dir.prepro,[subname,'_pre.set']);
    EEG = pop_loadset('filename',set_name);
    if strcmp(subname,'DMSS037')
        trialN = 738;
    else
        trialN = 720; % number of trials
    end
    %% time locked to delay onset

    if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
        markers = {'S 50'}; % delay onset
    else
        markers = {'50'}; % delay onset
    end

    timelimits = [-5.7 6];

    [sEEG,indices] = pop_epoch(EEG,markers,timelimits);
    rmvd = setdiff(1:trialN,indices);
    %% link behavoir, N =720 trials,Add new events to a loaded dataset 0.1 second after time-locking event
    trialID = 1:trialN;
    trialID(rmvd) = [];

    if sEEG.trials~=length(trialID)
        error('trial num mismatch')
    end

    nevents = length(sEEG.event);
    t = 0;
    for index = 1 : nevents
        if ischar(sEEG.event(index).type) && strcmpi(sEEG.event(index).type, markers)     % Add events relative to existing events
            t = t+1;
            markerStr = ['T' num2str(trialID(t))];
            sEEG.event(end+1) = sEEG.event(index); % Add event to end of event list
            % Specifying the event latency to be 0.1 sec before the referent event (in real data points)
            sEEG.event(end).latency = sEEG.event(index).latency + 0.1*sEEG.srate;
            sEEG.event(end).type = markerStr; % Make the type of the new event
        end
    end

    sEEG = eeg_checkset(sEEG, 'eventconsistency'); % Check all events for consistency

    %%

    new_name = [subname,'_epo.set'];
    pop_saveset(sEEG,'filename',new_name,'filepath',Dir.prepro);
end