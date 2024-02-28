clc
clear
% path to HMAX functions
base_dir = "/Users/kamp/PhD/spikevar/repo/SpikeVar";
cd(base_dir)
% add subfolders
addpath(genpath(base_dir));
% paths to stimuli
image_dir = "/Users/kamp/PhD/spikevar/stimuli_rn/";
% path to outputs
hmax_dir = fullfile(base_dir , "output", "hmax/");
vgg_dir = fullfile(base_dir , "output", "vgg16/");
neuro_dir = fullfile(base_dir, "output", "neuro/");
pls_dir = fullfile(base_dir, "output", "pls/");
behav_dir = fullfile(base_dir, "output", "behav/");
% data directories
data_dir = "/Users/kamp/PhD/spikevar/data/";
stim_dir = "/Users/kamp/PhD/spikevar/stimuli_rn";


%% Supplements
% Hippocampus
hpc_ids = find(pe_hpc_pls.lvlv_corrs.all_layers{1}>0);
plot_lvlv_corrs_all_layers(pe_hpc_pls)
plot_raincloud_stim_weights(pe_hpc_pls.stim_weights.late_layers, 'early', hpc_ids)
plot_raincloud_stim_weights(pe_hpc_pls.stim_weights.all_layers, 'all', hpc_ids)

% Amygdala
amyg_ids = find(pe_amyg_pls.lvlv_corrs.all_layers{1}>0);
plot_raincloud_stim_weights(pe_amyg_pls.stim_weights.early_layers, 'early', amyg_ids)
plot_raincloud_stim_weights(pe_amyg_pls.stim_weights.late_layers, 'late', amyg_ids)
plot_raincloud_stim_weights(pe_amyg_pls.stim_weights.all_layers, 'all', amyg_ids)

%%
% Compare neuron weights
load([data_dir + 'neuro/SpikeVar_VS_cells'])
load([data_dir + 'neuro/SpikeVar_MS_cells'])
load(hpc_sd_pls_file, 'all_neuron_ids');
compare_neuro_weights(pe_hpc_pls, sigCellListVS, sigCellListMS, hpc_ids, all_neuron_ids)