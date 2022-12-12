clear;

Dir.ana = 'C:\Users\Researcher\Documents\GitHub\WenWen\TimeEst\Ana_TE';
Dir.prepro = 'D:\rob\timeest\Prepro';
Dir.results = 'D:\rob\timeest\Results';
Dir.raw = 'D:\rob\timeest\Raw';
Dir.figs = 'D:\rob\timeest\FigsOutput';
Dir.toolbox = 'D:\intWM-E\toolbox';

subs = table;
tmp = dir(fullfile(Dir.raw,'*.vhdr'));
subs.name = cellfun(@(x)(['S' x(isstrprop(x,'digit'))]),{tmp.name}','UniformOutput',false);
subs.rawEEG = {tmp.name}';
subs.excluded = ismember(subs.name,{'S067','S069'});

%%
save(fullfile(Dir.ana,'subs.mat'),'subs','Dir')
