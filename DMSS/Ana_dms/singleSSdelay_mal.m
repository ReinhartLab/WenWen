function singleSSdelay_mal(sn,nIter,nfolds,IsOcci,IsBL2preDelay,IsCorretTrials,IsOverwrite)

load('subs.mat');
subname = subs.name{sn};

txtCell = {'','','';'_occi','_bl2preDelay','_corrTrials'};
outputFile = fullfile(Dir.results,[subname,'_ssDelay',txtCell{IsOcci+1,1},txtCell{IsBL2preDelay+1,2},txtCell{IsCorretTrials+1,3},'_mal.mat']);

if IsOverwrite==0 && isfile(outputFile)
    return
end

set_name = fullfile(Dir.prepro,[subname,'_clean.set']);
if isfile(set_name)
    EEG = pop_loadset('filename',set_name);

 %% load behavior

    csvFile = fullfile(Dir.beha,subs.csvFile{sn});
    M = readtable(csvFile);
    M = M(:,["block_num","ss_num","button_resp_rt","button_resp_corr"]);
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
        M = M(goodTrials,:);
    else
        error('trials mismatch')
    end
    %%

    if IsBL2preDelay
        EEG = pop_rmbase(EEG,[-200 0]);% in ms
    else
        bsWindow = {[-1.4 -1.2],[-2.6 -2.4],[-5 -4.8]};
        setsize = [1 2 4];
        for ss = 1:3
            tmpTrl = find(M.ss_num==setsize(ss));
            tmpEEG = pop_select(EEG,'trial',tmpTrl);

            tmpEEG = pop_rmbase(tmpEEG,bsWindow{ss}*1000);% pre trial baseline
            EEG.data(:,:,tmpTrl) = tmpEEG.data;
        end
    end

    EEG = pop_select(EEG,'time',[-1 4],'nochannel', {'TVEOG','LHEOG','RHEOG','BVEOG'});% [-0.2 4] relative to delay
    
    if IsOcci
        occiChans = {'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2'};
        EEG = pop_select(EEG,'channel',occiChans);
    end
   

    %% decode
    clear dcd

    if IsCorretTrials ==1
        tmpID =  M.button_resp_corr ==1;
    else
        tmpID = ones(height(M),1);% all trials
    end

    EEG = pop_select(EEG,'trial',find(tmpID));
    ssArray = M.ss_num (find(tmpID),:);
    ssArray(ssArray==4)=3;

    dat = permute(EEG.data,[3 1 2]);% trials*chan*time;

    nSamp = size(dat,3);

    dec_acc = nan(nIter,nSamp);
    dat = smoothdata(dat,3,'gaussian',40);

    parfor iter = 1:nIter
        dec_acc(iter,:) = mahal_nominal_kfold_b(dat,ssArray,nfolds);
    end

    dcd.ACC = mean(dec_acc,'omitnan');

    dcd.Dtime = [EEG.xmin EEG.xmax];
    dcd.times = EEG.times;
    dcd.nSamps = nSamp;
    %%
    save(outputFile,'dcd','-v7.3')
end
