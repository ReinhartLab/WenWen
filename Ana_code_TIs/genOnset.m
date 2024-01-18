clear;
load('subsinfo.mat')

%%

for sn = 1:subN
    subname = subs.name{sn};
    %--------speficify subjects
    % subname = {'S01'};

    %% FFA localizer
    curTask = 'FFA';
    condStr = {'Image','Scrambled'};
    mat_name = dir(fullfile(Dir.beha.(curTask),subname,[subname '*',curTask,'.mat'])).name;
    ffa = load(fullfile(Dir.beha.(curTask),subname,mat_name));

    % check timing accuracy
    clear tmp
    tmp(1) = (ffa.trials.curTR_s(end)-1)+(ffa.params.stimuli+ffa.params.response)/ffa.params.TR+ffa.trials.itiTR(end)+ffa.params.overrun(2)/ffa.params.TR;
    tmp(2) = sum(diff(ffa.trials.curTR_s) == ffa.trials.itiTR(1:end-1)+(ffa.params.stimuli+ffa.params.response)/ffa.params.TR + 1);
    if tmp(1) == ffa.ideal_Vol && tmp(2) == ffa.params.trlNperRun-1
        outMat =table;
        tmfun = @(x1,x2)[x1, '_loc',num2str(x2)];
        outMat.trial_type = cellfun(tmfun,condStr(ffa.trials.scrambled+1)',num2cell(ffa.trials.faceposi),'UniformOutput',false);

        outMat.onset = (ffa.trials.curTR_s-1)*ffa.params.TR-ffa.params.preScan;%onset time in s
        outMat.duration = ones(height(ffa.trials),1)*ffa.params.stimuli;%duration
        outMat.weight = ones(height(ffa.trials),1); % parametric
        destFile = fullfile(Dir.timing,curTask,[subname '.csv']);
        writetable(outMat,destFile);

    else
        error(sprintf('%s timing error',curTask))
    end

    %% MT localizer
    curTask = 'MT';
    condStr = {'Motion','Static'};
    mat_name = dir(fullfile(Dir.beha.(curTask),subname,[subname '*',curTask,'.mat'])).name;
    mt = load(fullfile(Dir.beha.(curTask),subname,mat_name));

    onset = (mt.trials.curTR_s-1)*mt.params.TR-mt.params.preScan;%in s

    % check timing accuracy
    clear tmp
    tmp(1) = (mt.trials.curTR_s(end)-1)+(mt.params.stimuli+mt.params.response)/mt.params.TR+mt.trials.itiTR(end)+mt.params.overrun(2)/mt.params.TR;
    tmp(2) = sum(diff(mt.trials.curTR_s) == mt.trials.itiTR(1:end-1)+(mt.params.stimuli+mt.params.response)/mt.params.TR + 1);
    if tmp(1) == mt.ideal_Vol && tmp(2) == mt.params.trlNperRun-1
        outMat =table;
        tmfun = @(x1,x2)[x1, '_loc',num2str(x2)];
        outMat.trial_type = cellfun(tmfun,condStr(mt.trials.static+1)',num2cell(mt.trials.motionposi),'UniformOutput',false);

        outMat.onset = (mt.trials.curTR_s-1)*mt.params.TR-mt.params.preScan;%onset time in s
        outMat.duration = ones(height(mt.trials),1)*mt.params.stimuli;%duration
        outMat.weight = ones(height(mt.trials),1); % parametric
        destFile = fullfile(Dir.timing,curTask,[subname '.csv']);
        writetable(outMat,destFile);

    else
        error(sprintf('%s timing error',curTask))
    end

    %% OBA localizer
    curTask = 'OBA';
    condStr = {'Image','Scrambled'};
    mat_name = dir(fullfile(Dir.beha.(curTask),subname,[subname '*',curTask,'.mat'])).name;
    oba = load(fullfile(Dir.beha.(curTask),subname,mat_name));

    onset = (oba.trials.curTR_s-1)*oba.params.TR-oba.params.preScan;%in s

    % check timing accuracy
    clear tmp
    tmp(1) = (oba.trials.curTR_s(end)-1)+(oba.params.stimuli+oba.params.response)/oba.params.TR+oba.trials.itiTR(end)+oba.params.overrun(2)/oba.params.TR;
    tmp(2) = sum(diff(oba.trials.curTR_s) == oba.trials.itiTR(1:end-1)+(oba.params.stimuli+oba.params.response)/oba.params.TR + 1);
    if tmp(1) == oba.ideal_Vol && tmp(2) == oba.params.trlNperRun-1
        outMat =table;
        tmfun = @(x1,x2)[x1, '_loc',num2str(x2)];
        outMat.trial_type = cellfun(tmfun,condStr(oba.trials.scrambled+1)',num2cell(oba.trials.vehiposi),'UniformOutput',false);

        outMat.onset = (oba.trials.curTR_s-1)*oba.params.TR-oba.params.preScan;%onset time in s
        outMat.duration = ones(height(oba.trials),1)*oba.params.stimuli;%duration
        outMat.weight = ones(height(oba.trials),1); % parametric
        destFile = fullfile(Dir.timing,curTask,[subname '.csv']);
        writetable(outMat,destFile);

    else
        error(sprintf('%s timing error',curTask))
    end

    %% IBL task

    % curTask = 'fMRI';
    %
    % load(fullfile(Dir.beha.IBL,subname,[subname,'.mat']))
    %     fullsubName = subs{contains(subs,subname,'IgnoreCase',true)};
    %     subname = fullsubName(1:end-5);
    %     load(fullfile(curTaskFolder,[subname,'_fmri.mat']));
    %
    %     condStr = {'predicted','random'};
    %
    %     %--------------block design
    %
    %     for run_i = 1:8
    %         trials = stimA(run_i).all;
    %         onset = ([trials.curTR_s]'-1-params.preScan)*params.TR;%in s
    %
    %         blockS = 1:params.trialsPerblock:params.trialsPerRun;
    %         blockE = [1:4]*params.trialsPerblock;
    %
    %         clear outMat
    %         outMat(:,1) = onset(blockS);%onset time
    %         outMat(:,2) = onset(blockE)+[trials(blockE).itiTR]'*params.TR - onset(blockS)+1*params.TR;%duration
    %         outMat(:,3) = ones(4,1); % parametric
    %
    %         for cond_i = 1:2
    %             niiFolder = fullfile(Dir.data,subname,['day',num2str(day_i)],['FBE',num2str(run_i)]);
    %             destFile = fullfile(niiFolder,[condStr{cond_i},'_blocks.txt']);
    %             tmp = contains({trials(blockS).cond},condStr{cond_i});
    %             tmpOut= outMat(tmp,:);
    %             save(destFile,'tmpOut','-ascii');
    %         end
    %     end
    %
    %     % ---------------event-analysis: orientation
    %
    %     for run_i = 1:8
    %         trials = stimA(run_i).all;
    %         onset = ([trials.curTR_s]'-1-params.preScan)*params.TR;%in s
    %         niiFolder = fullfile(Dir.data,subname,['day',num2str(day_i)],[curTask,num2str(run_i)]);
    %
    %         tmpLocs = cat(1,trials.locs);
    %
    %         % to check timing accuracy
    %         if double(trials(end).curTR_s+trials(end).itiTR == 201)% postblock+overrun = 9TR
    %             for loci = 1:4
    %                 tmp_ori_array = tmpLocs(:,loci).*[trials.bin1]'+ ~logical(tmpLocs(:,loci)).*[trials.bin2]';
    %                 for ori = 1:params.OriN
    %                     clear outMat
    %
    %                     outMat(:,1) = onset(tmp_ori_array==ori);%onset time
    %                     outMat(:,2) = ones(sum(tmp_ori_array==ori),1)*params.TR;%duration
    %                     outMat(:,3) = ones(sum(tmp_ori_array==ori),1); % parametric
    %
    %                     destFile = fullfile(niiFolder,['FBE',num2str(run_i),'_loc',num2str(loci),'_ori', num2str(ori) '.txt']);
    %                     save(destFile,'outMat','-ascii');
    %                 end
    %             end
    %         else
    %             error('timing error')
    %         end
    %     end
    %
    %     % ---------------event-analysis: event onset
    %
    %     for run_i = 1:8
    %         trials = stimA(run_i).all;
    %         onset = ([trials.curTR_s]'-1-params.preScan)*params.TR;%in s
    %         niiFolder = fullfile(Dir.data,subname,['day',num2str(day_i)],[curTask,num2str(run_i)]);
    %
    %         clear outMat
    %         outMat(:,1) = onset;%onset time
    %         outMat(:,2) = ones(params.trialsPerRun,1)*params.TR;%duration
    %         outMat(:,3) = ones(params.trialsPerRun,1); % parametric
    %
    %         for cond_i = 1:2
    %             destFile = fullfile(niiFolder,[condStr{cond_i},'_ERonset.txt']);
    %             tmp = contains({trials.cond},condStr{cond_i});
    %             tmpOut= outMat(tmp,:);
    %             save(destFile,'tmpOut','-ascii');
    %         end
    %     end
    %
    %
    %     % ---------------event-analysis: parametric for pred, 1 for rand
    %
    %     for run_i = 1:6
    %         trials = stimA(run_i).all;
    %         onset = ([trials.curTR_s]'-1-params.preScan)*params.TR;%in s
    %         niiFolder = fullfile(Dir.data,subname,['day',num2str(day_i)],[curTask,num2str(run_i)]);
    %
    %         clear outMat
    %         outMat(:,1) = onset;%onset time
    %         outMat(:,2) = ones(params.trialsPerRun,1)*params.TR;%duration
    %         outMat(:,3) = ones(params.trialsPerRun,1); % parametric
    %
    %         for cond_i = 1:2
    %             destFile = fullfile(niiFolder,[condStr{cond_i},'_ERonset_para.txt']);
    %             tmp = contains({trials.cond},condStr{cond_i});
    %             tmpOut= outMat(tmp,:);
    %             if cond_i == 1
    %                 tmpOut(:,3) = repmat([1:-1/16:1/16]',size(tmpOut,1)/16,1);
    %             end
    %             save(destFile,'tmpOut','-ascii');
    %         end
    %     end

end