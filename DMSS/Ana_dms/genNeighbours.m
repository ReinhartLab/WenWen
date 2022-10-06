clear;
addpath('D:\Toolbox\fieldtrip-20211209')
ft_defaults
%%
cfg = [];
elec = ft_read_sens('standard_1020.elc');
ft_plot_sens(elec, 'label', 'yes', 'elecshape', 'disc', 'elecsize', 10, 'facecolor', [0.8 0.8 1.0])
camlight headlight

%%
cfg = [];
cfg.elec = elec;
cfg.method = 'distance';
cfg.feedback = 'yes';
cfg.neighbourdist = 50;
neighbours = ft_prepare_neighbours(cfg);

save neighbours neighbours