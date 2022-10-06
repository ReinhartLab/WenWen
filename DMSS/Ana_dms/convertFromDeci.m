
Deci.Folder.Raw = 'D:\dms_fb\processed\Artifact';

Deci.SubjectList = cellfun(@(c) strsplit(c,'.'),CleanDir(Deci.Folder.Raw),'un',0);
Deci.SubjectList = unique(cellfun(@(c) c{1},Deci.SubjectList,'un',0));

fakeUI = figure;
fakeUI.UserData = Deci.SubjectList;
fakeUI.Visible =  'off';
dc_select_labels(fakeUI,[],Deci.SubjectList);
waitfor(findall(0,'Name','Select Labels'),'BeingDeleted','on');
Deci.SubjectList = fakeUI.UserData;
close(fakeUI);

Deci.Plot.CondTitle     = {'Set Size 1 Delay' 'Set Size 2 Delay' 'Set Size 4 Delay' 'SS1 Correct Rsp' 'SS2 Correct Rsp' 'SS4 Correct Rsp'};

Chois = CleanDir('D:\dms_fb\processed_Delay1\Analysis\Volt_Raw\dmss01\Delay Onset\Set Size 1 Delay');
info.variable = 'raw';

Deci.Plot.Lock = {'First Stim Onset' 'Delay Onset' 'Rsp Onset'}; 

Deci.Folder.Analysis = 'D:\dms_fb\processed_Delay1\Analysis';
%%

for  subject_list = 22:length(Deci.SubjectList)

    display(['Loading Plottor for Subject #' num2str(subject_list) ': ' Deci.SubjectList{subject_list}]);
    clear Subjects

    for Lock = 1:length(Deci.Plot.Lock)

        for Conditions = 1:length(Deci.Plot.CondTitle)
            for Channel = 1:length(Chois)

                load([Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock{Lock} filesep Deci.Plot.CondTitle{Conditions} filesep Chois{Channel}],'raw');



                Chans{Channel} = raw;
                Chans{Channel}.label = Chois(Channel);

                clear raw;
            end

            if isfield(Chans{Channel},'trllength')
                trllen(subject_list,Conditions) = Chans{Channel}.trllength;
            else
                trllen(subject_list,Conditions) = nan;
            end

            acfg.parameter = 'avg';
            acfg.appenddim = 'chan';
            Subjects{Conditions} = rmfield(ft_appendtimelock(acfg,Chans{:}),'cfg');
            Subjects{Conditions}.dimord = 'chan_time';
            Subjects{Conditions}.cond = Deci.Plot.CondTitle{Conditions};

            Subjects{Conditions}.label = cellfun(@(x)x(1:end-4),Subjects{Conditions}.label,'un',0);

        end

        subDir = [Deci.Folder.Analysis filesep 'Volt_Raw' filesep Deci.SubjectList{subject_list}  filesep Deci.Plot.Lock{Lock} ];
        tmp_name = 'appended.mat';
        disp(subDir)
        save(fullfile(subDir,tmp_name),'Subjects','-v7.3')

    end

    clear Chans;
end