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

%% 1. Neuro: Compute spike permutation entropy
% 1.1. Prepare the binned spiking data from the raw data
spike_file = [neuro_dir + "SpikeVar_spike_data.mat"];
if isfile(spike_file)
    fprintf("Using existing binned spiking data from output/neuro/SpikeVar_spike_data.mat.\n")
    load(spike_file, 'spike_data_table');
else 
    fprintf("Extracting binned spiking data.\n")
    spike_data_table = a_prepare_spike_data(base_dir, data_dir, stim_dir);   
end

% 1.2. Calculate spike permutation entropy
%spike_data_table = b_compute_spike_pe(spike_data_table, neuro_dir); 

% alternatively load precomputed PE table
load([neuro_dir + 'SpikeVar_spike_data_PE.mat'], 'spike_data_table')

% 1.3. compute spike rate, average PE etc.
spike_summary = c_get_spike_summay(spike_data_table, neuro_dir);

%% 2. Stimuli: Compute feature estimates using HMAX and VGG
% load behavioral information and stimuli information
load([data_dir + 'behavior/SpikeVar_learn_recog_info.mat'], 'fulltable_learn', 'fulltable_recog');
stim_info = create_stim_info(image_dir, fullfile(base_dir, "output/")); 

% 2.1. Get HMAX estimates
hmax_file = [hmax_dir + 'HMAX_estimates.mat']; 
if isfile(hmax_file)
    fprintf("Using precomputed HMAX activation estimates.\n")
    load(hmax_file, 'HMAX_estimates');
else 
    fprintf("Estimating HMAX feature activations.\n")
    fprinft("This may take a couple of minutes.\n")
    HMAX_estimates = a_get_hmax_estimates(image_dir, hmax_dir);
end 

% 2.2. Run pca HMAX features estimates 
fprintf("Running PCA on HMAX feature estimates.\n")
hmax_pc1_table = b_hmax_pca(HMAX_estimates, stim_info, fulltable_learn, fulltable_recog, hmax_dir);
%load('/Users/kamp/PhD/spikevar/repo/SpikeVar/stimuli/old_version/original_results_renamed/renamed_hmax.mat', 'hmax_pc1_table')

% 2.3. Get VGG16 estimates
% We're using the precomputed VGG16 estimates. To compute them again, run
% stimuli/c_get_vgg_estimates.py externally

% 2.4. Run pca on VGG16 estimates
fprintf("Running PCA on VGG16 feature estimates.\n")
vgg_dir_or = fullfile(base_dir , "stimuli", "old_version", "original_results_renamed/");
vgg_pc1_table = d_vgg16_pca(vgg_dir_or, stim_info, fulltable_learn, fulltable_recog);


%% 3. Partial Least Squares 
% 3.1. Run PLS analysis on hippocampus and amygdala data
pe_hpc_pls_file = [pls_dir + 'SpikeVar_PE_PLS_hippocampus_results.mat'];
pe_amyg_pls_file = [pls_dir + 'SpikeVar_PE_PLS_amygdala_results.mat'];
spikenum_hpc_pls_file = [pls_dir + 'SpikeVar_SpikeNum_PLS_hippocampus_results.mat'];

if 1==0 %isfile(pe_hpc_pls_file)
    fprintf("Using precomputed PLS estimates.\n")
    load(pe_hpc_pls_file, 'pls_results'); pe_hpc_pls = pls_results;
    load(pe_amyg_pls_file, 'pls_results'); pe_amyg_pls = pls_results;
    load(spikenum_hpc_pls_file, 'pls_results'); spikenum_hpc_pls = pls_results;
else
    % compute pls based on permutation entropy
    pe_hpc_pls = a_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, 'hippocampus', 'PE');
    pe_amyg_pls = a_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, 'amygdala', 'PE');
    % compute pls based on the spike number
    spikenum_hpc_pls = a_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, 'hippocampus', 'SpikeNum');
    % spikenum_amyg_pls = a_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, 'amygdala', 'SpikeNum');
end

% 3.2. Analyse PLS results 
% 3.2.1. Plot latent correlations of the Hippocampus
hpc_ids = find(pe_hpc_pls.lvlv_corrs.all_layers{1}>0);
plot_lvlv_corrs(pe_hpc_pls.lvlv_corrs.late_layers{1}, pe_hpc_pls.lvlv_corr_lims.late_layers{1}, ... 
    hpc_ids, 'Hippocampus PE vs late layers')

% 3.2.2. Plot stimulus weights of the Hippocampus
plot_stim_weights(pe_hpc_pls.stim_weights.late_layers, hpc_ids)
plot_raincloud_stim_weights(pe_hpc_pls.stim_weights.late_layers, 'late', hpc_ids)

% 3.2.4. save absolute stim weights to compute repeated measure models later
save_abs_stim_weights(pe_hpc_pls.stim_weights, hpc_ids, pls_dir)

% 3.3. run visual stimulus weight mixing
% b_vis_stim_weight_mixing(spike_data_table, hmax_pc1_table, vgg_pc1_table)

%% 4. Behavioral Models
% 4.1. Prepare behavioral data
load([data_dir + 'behavior/SpikeVar_learn_recog_info.mat'], 'fulltable_learn', 'fulltable_recog');
pat_ids = unique(spike_data_table.Participant);
behav_data = a_extract_behav_info(fulltable_learn, fulltable_recog, pat_ids, behav_dir);

% 4.2. Convert to ranks
load(pe_hpc_pls_file, 'trial_nums', 'neuron_nums');
rank_table = b_convert2ranks(pe_hpc_pls, pe_amyg_pls, spikenum_hpc_pls, ...
    behav_data, neuron_nums, trial_nums, spike_summary);

% 4.3. Run pca on performance data
rel_inds.hpc = find(pe_hpc_pls.lvlv_corrs.all_layers{1}~=0);
rel_inds.amyg = find(pe_amyg_pls.lvlv_corrs.late_layers{1}(rel_inds.hpc)~=0);
rank_table = c_behav_pca(rank_table, rel_inds);
save('rank_table.mat', 'rank_table')

% 4.4. Run linear model
d_lm_layer_vs_perf(rank_table)

% 4.5. Control analysis 
e_control_analysis(rank_table)






