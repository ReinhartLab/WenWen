function loadRaw(sn)

%%
load('subs.mat')
if subs.rawEEG(sn)

    subname = subs.name{sn};

    subDir = fullfile(Dir.raw,subname);

    files = dir(subDir);
    tmp_idx = find(contains({files.name},'.vhdr'));
    eeg_name = {files(tmp_idx).name};

    if length(eeg_name) == 1
        EEG = pop_loadbv(subDir,eeg_name{1});
    else
        [ALLEEG,~,~] = eeglab;

        for i = 1:length(eeg_name)
            EEG = pop_loadbv(subDir,eeg_name{i});
            [ALLEEG,~] = eeg_store(ALLEEG,EEG);
        end

        EEG = pop_mergeset(ALLEEG,1:length(eeg_name));
    end

    set_name = [subname,'_raw.set'];
    pop_saveset(EEG,'filename',set_name,'filepath',Dir.prepro);

    %%
    EEG = pop_loadset('filename',set_name,'filepath',Dir.prepro);
    EEG = pop_resample(EEG,500);
    EEG = pop_eegfiltnew(EEG,0.5,40);
    EEG = pop_chanedit(EEG, 'lookup',fullfile('D:\intWM-E\toolbox\eeglab2021.1\plugins\dipfit4.3\standard_BESA','standard-10-5-cap385.elp'));

    % online eeg recording use TP10 for reference, we add a dummy TP10 and re-reference to the average of TP9&10
    chanID = 64; % the original ch64 is stimtrak
    EEG.nbchan      = chanID;
    EEG.data(chanID,:)  = zeros(1,EEG.pnts);
    EEG.chanlocs(chanID).label = 'TP10';
    EEG = pop_chanedit(EEG, 'changefield',{chanID 'labels' 'TP10'},'lookup',fullfile('D:\intWM-E\toolbox\eeglab2021.1\plugins\dipfit4.3\standard_BESA','standard-10-5-cap385.elp'), ...
        'eval','chans = pop_chancenter( chans, [],[]);');

    EEG = pop_reref( EEG, {'TP9','TP10'});

    set_name = [subname,'_rawFiltered.set'];
    pop_saveset(EEG,'filename',set_name,'filepath',Dir.prepro);

end