clear
close all
load('subs.mat');

parfor sn = 1:height(subs)
    if subs.rawEEG(sn)
        subname = subs.name{sn};
        set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
        EEG = pop_loadset('filename',set_name);

        %%
        %         if sum(strcmp({EEG.event.type},'S 50'))>0
        if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
            stim1 = {'S 11','S 21', 'S 31'};
        else
            stim1 = {'11','21','31'};
        end
      
        bEEG = pop_select(EEG,'time',[-5.2 2]);% re-epoch into shorter segments
        timelimits = [-0.2 0];% pre-stimuli baseline
        [bEEG,indices_b] = pop_epoch(bEEG,stim1,timelimits);

        if EEG.trials~=bEEG.trials
            [EEG.trials bEEG.trials]
            error('trials missing')
        end

        %% baseline subtraction
        EEG.data = bsxfun(@minus,EEG.data,mean(bEEG.data,2));

        set_name = fullfile(Dir.prepro,[subname,'_delay.set']);
        pop_saveset(EEG,'filename',set_name);

    end
end