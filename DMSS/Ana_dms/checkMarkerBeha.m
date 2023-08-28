function checkMarkerBeha(EEG,M)

%%
behaCode = cell(height(M),5);
if sum(strcmp({EEG.event.type},'S 50'))>=sum(strcmp({EEG.event.type},'50'))
    stim = 1;%{'S 11','S 21', 'S 31'};
else
    stim = 2;%{'11','21','31'};
end

for trl_i = 1:height(M)
    if table2array(M(trl_i,"ss_num"))==1
        ssCode = {'11','12','52'};
    elseif table2array(M(trl_i,"ss_num"))==2
        ssCode = {'25','26','53'};
    else
        ssCode = {'43','44','54'};
    end

    if table2array(M(trl_i,"button_resp_corr"))==1
        respCode = {'64'};
    else
        respCode = {'65'};
    end

    if table2array(M(trl_i,"button_resp_keys"))==1
        keyCode = {'61'};
    else
        keyCode = {'62'};
    end
    curMarkers = [ssCode,respCode,keyCode];
    if stim == 1
        curMarkers = cellfun(@(x)['S' sprintf('%3s',num2str(x))],curMarkers,'UniformOutput',false);
    end
    behaCode(trl_i,:) = curMarkers(:)';
end

for trl_i = 1:EEG.trials
    if sum(ismember(behaCode(trl_i,:),EEG.epoch(trl_i).eventtype))<3
        error('mismatch triggers at trial %d',trl_i)
    end
end
