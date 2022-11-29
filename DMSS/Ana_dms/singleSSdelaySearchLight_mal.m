function singleSSdelaySearchLight_mal(sn,nIter,nfolds,IsBL2preDelay,IsCorretTrials,IsOverwrite)

load('subs.mat');
load neighbours.mat
subname = subs.name{sn};

txtCell = {'','';'_bl2preDelay','_corrTrials'};
outputFile = fullfile(Dir.results,[subname,'_ssDelayChan',txtCell{IsBL2preDelay+1,1},txtCell{IsCorretTrials+1,2},'_mal.mat']);

if isfile(outputFile) && IsOverwrite ==0
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

    %%

    neighbours = neighbours(ismember({neighbours.label},{EEG.chanlocs.labels}));

    for c = 1:length({neighbours.label})
        curChan = neighbours(c).label;
        curNB = neighbours(c).neighblabel;
        [neighbours(c).ismemb,neighbours(c).neighbID]= ismember(curNB,{EEG.chanlocs.labels});
        neighbours(c).neighbID = neighbours(c).neighbID(neighbours(c).ismemb);
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
    nChan = length({neighbours.label});
    dec_acc = nan(nIter,nChan,nSamp);

    for c = 1:nChan
        cdata = dat(:,neighbours(c).neighbID,:);% trials*chan*time;
        parfor iter = 1:nIter
            dec_acc(iter,c,:) = mahal_nominal_kfold_b(cdata,ssArray,nfolds);
        end
    end

    dcd.ACC = mean(dec_acc,'omitnan');

    dcd.Dtime = [EEG.xmin EEG.xmax];
    dcd.times = EEG.times;
    dcd.nSamps = nSamp;
    dcd.neighbours = neighbours;
    %%
    save(outputFile,'dcd','-v7.3')

end
