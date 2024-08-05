clc
clear
% path to HMAX functions
base_dir = "XXXX"; % enter the directory of the repository here
cd(base_dir)
% add subfolders
addpath(genpath(base_dir));

% path to outputs
hmax_dir = fullfile(base_dir , "output", "hmax/");
vgg_dir = fullfile(base_dir , "output", "vgg16/");
neuro_dir = fullfile(base_dir, "output", "neuro/");
pls_dir = fullfile(base_dir, "output", "pls/");
behav_dir = fullfile(base_dir, "output", "behav/");
perm_dir = fullfile(base_dir, "output", "perm/");

% data directories
data_dir = fullfile(base_dir , "data/");
stim_dir = fullfile(base_dir , "stimuli/");

% paths to stimuli
image_dir = fullfile(data_dir , "stimuli_rn");
%% 1. Neuro: Compute spike permutation entropy
% 1.1. Prepare/load the binned spiking data from the raw data
spike_file = [neuro_dir + "SpikeVar_spike_data_PE_or.mat"];
load(spike_file, 'spike_data_table')

% 1.2. compute spike rate, average PE etc.
spike_summary = c_get_spike_summay(spike_data_table, neuro_dir);

%% 2. Stimuli: Compute feature estimates using HMAX and VGG
% load behavioral information and stimuli information
load([data_dir + 'SpikeVar_learn_recog_info.mat'], 'fulltable_learn', 'fulltable_recog');
stim_info = create_stim_info(image_dir, fullfile(base_dir, "output/")); 

% 2.1. Get HMAX estimates
hmax_file = [hmax_dir + 'HMAX_estimates_or.mat']; 
fprintf("Using precomputed HMAX activation estimates.\n")
load(hmax_file, 'HMAX_estimates');

% 2.2. Run pca HMAX features estimates 
fprintf("Running PCA on HMAX feature estimates.\n")
hmax_pc1_table = b_hmax_pca(HMAX_estimates, stim_info, fulltable_learn, fulltable_recog, hmax_dir);
%load('/Users/kamp/PhD/spikevar/repo/SpikeVar/stimuli/old_version/original_results_renamed/renamed_hmax.mat', 'hmax_pc1_table')

% 2.3. Get VGG16 estimates
% We're using the precomputed VGG16 estimates. To compute them again, run
% stimuli/c_get_vgg_estimates.py externally

% 2.4. Run pca on VGG16 estimates
fprintf("Running PCA on VGG16 feature estimates.\n")
vgg_dir_or = fullfile(vgg_dir , "original_results/");
vgg_pc1_table = d_vgg16_pca(vgg_dir_or, stim_info, fulltable_learn, fulltable_recog);

%% 3. Partial Least Squares 
% compute pls 
% repeat for multiple permutations
[pls_results, permutations, info] = perm_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table);
save([perm_dir + 'pls_permutation_results.mat'], 'pls_results', 'info', 'permutations')
