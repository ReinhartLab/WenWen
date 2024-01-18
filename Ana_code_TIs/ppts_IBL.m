clear;
Dir.ana = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code';
Dir.beha.FFA = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Behav/FFA';
Dir.beha.OBA = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Behav/OBA';
Dir.beha.MT = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Behav/MT';
Dir.beha.IBL = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Behav/fMRI';
Dir.timing = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Ana_code/timing';
Dir.figs = '/projectnb/viscog01/Wen/IBL_TI_fMRI/Figs';
Dir.bids = '/projectnb2/viscog01/Wen/IBL_TI_fMRI/BIDS/';
% if ~isfolder(Dir.results),mkdir(Dir.results),end
if ~isfolder(Dir.timing),mkdir(Dir.timing),end

%% using bids dir to obtain subject names 

files = dir(fullfile(Dir.bids,'sub-*'))
names = {files.name};

names = cellfun(@(x)extractAfter(x,numel(x)-3),names,'UniformOutput',false);
fid = fopen(fullfile(Dir.bids,'subj_all.txt'), 'w');
fprintf(fid, '%s\n', names{:});
fclose(fid);

%%
subN = length(names)
subs = table;
subs.name = names';

excludeSubj = {'S60'};

subs.exclude = contains(subs.name,excludeSubj);

save(fullfile(Dir.ana,'subsinfo.mat'),'subN','Dir','subs')
