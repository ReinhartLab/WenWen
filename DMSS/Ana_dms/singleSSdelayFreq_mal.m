function singleSSdelayFreq_mal(sn,nIter,nfolds,IsOcci,IsdePhase,IsBL2preDelay,IsCorretTrials,IsOverwrite)

load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','','';'_occi','_dephase','_bl2preDelay','_corrTrials'};
outputFile =fullfile(Dir.results,[subname,'_ssDelayFreq',txtCell{IsOcci+1,1},txtCell{IsdePhase+1,2},txtCell{IsBL2preDelay+1,3},txtCell{IsCorretTrials+1,4},'_mal.mat']);

if IsOverwrite==0 && isfile(outputFile)
    return
end


set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
if isfile(set_name)
    EEG = pop_loadset('filename',set_name);% baseline corrected to pre-stim
    %% load behavior

    csvFile = fullfile(Dir.beha,subs.csvFile{sn});
    M = readtable(csvFile);
    M = M(:,["block_num","ss_num","type","button_resp_rt","button_resp_corr"]);
    M(1:end-1,{'button_resp_corr','button_resp_rt'}) = M(2:end,{'button_resp_corr','button_resp_rt'});

    M(isnan(M.ss_num),:) = [];
    if sum(isnan(M.button_resp_rt)) >0
        M(isnan(M.button_resp_rt),["button_resp_rt","button_resp_corr"]) = array2table([0 0]);% no response = wrong response
    end

    tmp = find(contains({EEG.event.type},'T'));
    if length(tmp) ==EEG.trials
        goodTrials = {EEG.event(tmp).type};
        goodTrials = cellfun(@(x)(str2num(x(2:end))),goodTrials,'UniformOutput',false);
        goodTrials = cell2mat(goodTrials);
    end
    M = M(goodTrials,:);
    %%
    EEG = pop_select(EEG,'nochannel', {'TVEOG','LHEOG','RHEOG','BVEOG'});
    if IsBL2preDelay
        EEG = pop_rmbase(EEG,[-400 0]);% in ms
    else
        bsWindow = {[-1.6 -1.2],[-2.8 -2.4],[-5.2 -4.8]};
        setsize = [1 2 4];
        for ss = 1:3
            tmpTrl = find(M.ss_num==setsize(ss));
            tmpEEG = pop_select(EEG,'trial',tmpTrl);

            tmpEEG = pop_rmbase(tmpEEG,bsWindow{ss}*1000);% pre trial baseline
            EEG.data(:,:,tmpTrl) = tmpEEG.data;
        end
    end
    if IsOcci
        occiChans = {'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2'};
        EEG = pop_select(EEG,'channel',occiChans);
    end

    %% decode
    clear dcd

    timeIdx = EEG.times<=4000 & EEG.times>=-1200;

    dcd.freqs = [[1:2:40];[1:2:40]+2];
    dcd.freqN = length(dcd.freqs);
    dcd.srate = 500;

    if IsCorretTrials ==1
        tmpID =  M.button_resp_corr ==1;
    else
        tmpID = ones(height(M),1);% all trials
    end

    EEG = pop_select(EEG,'trial',find(tmpID));
    ssArray = M.ss_num (find(tmpID),:);
    ssArray(ssArray==4)=3;

    dat = EEG.data;% chan*time*trial

    if IsdePhase
        dat = bsxfun(@minus,dat,mean(dat,3)); % remove ERP
    end

    nSamp = size(dat,2);

    dec_acc = nan(nIter,nSamp);

    for fi = 1:dcd.freqN
        dat_dec = dat;% chan*time*trial

        parfor c = 1:size(dat_dec,1)
            dat_dec(c,:,:) = abs(hilbert(eegfilt(squeeze(dat_dec(c,:,:))',dcd.srate,dcd.freqs(1,fi),dcd.freqs(2,fi))')).^2; % hilbert: time must be in the first dimension
        end
        dat_dec = permute(dat_dec,[3 1 2]);% trials*chan*time;

        parfor iter = 1:nIter
            dec_acc(iter,:) = mahal_nominal_kfold_b(dat_dec,ssArray,nfolds);
        end

        dcd.ACC{fi} = mean(dec_acc,'omitnan');
        dcd.ACC{fi} = dcd.ACC{fi}(timeIdx);
    end


    dcd.times = EEG.times(timeIdx);
    dcd.Dtime = [dcd.times(1) dcd.times(end)];

    dcd.nSamps = length(dcd.times);
    %%
    save(outputFile,'dcd','-v7.3')

end
