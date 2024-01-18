ppts_IBL % set up directory


%% individual
for sub_i = 1:subN

    subname = subs.name{sub_i};
    
    try
        PlotSingleBeha(subname)
        %     PlotSingleBeha('S02')
    catch
        fprintf('Skipped %s\n',subname)
    end

end
