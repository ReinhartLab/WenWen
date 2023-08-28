clear;

% Dir.ana = 'E:\dmss_OTVV\Ana_dms';
% Dir.beha = 'E:\dmss_OTVV\dmss behv csv files';
% Dir.prepro = 'E:\dmss_OTVV\prepro';
% Dir.results = 'E:\dmss_OTVV\results';
% Dir.raw = 'E:\dmss_OTVV\raw';

Dir.ana = 'C:\Users\Peyton\Documents\GitHub\WenWen\DMSS\Ana_dms';
Dir.beha = 'D:\dmss_OTVV\dmss behv csv files';
Dir.prepro = 'D:\dmss_OTVV\prepro';
Dir.results = 'D:\dmss_OTVV\results';
Dir.raw = 'D:\dmss_OTVV\raw';
Dir.figs = 'D:\dmss_OTVV\FigsOutput';

subs = table;
tmp = dir(Dir.beha);
subs.name = {tmp([tmp.isdir]==0).name}';
subs.csvFile = subs.name;

dmsfun = @(x,b)(x(1:b-1));
subs.name = cellfun(@(x,b) dmsfun(x,b),subs.name,strfind(subs.name,'_dms_ss_vv_v8'),'UniformOutput',false);
%%
rawsubs = {dir(Dir.raw).name};
rawsubs = rawsubs(contains(rawsubs,'dms','IgnoreCase',true));
subs.rawEEG = contains(subs.name,rawsubs,'IgnoreCase',true);

%%
groupStr = {'Young','Old'};
subs.group = double(contains(subs.name,'DMSS'))+1;
subs.groupStr = groupStr([subs.group])';

%%

rmvSubs = {'04','23','28'};
tmp_idx = ismember(subs.name,rmvSubs);
subs(tmp_idx,:) = [];

excludeSubj = {'dmss24','dmss09'};
subs.exclude = ismember(subs.name,excludeSubj);

%%
save(fullfile(Dir.ana,'subs.mat'),'subs','Dir')
