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


% data directories
data_dir = "/Users/kamp/PhD/spikevar/data/";
stim_dir = "/Users/kamp/PhD/spikevar/stimuli_rn";

%% Neuro: Compute spike permutation entropy
% 1. Prepare the binned spiking data from the raw data
spike_file = [neuro_dir + "SpikeVar_spike_data.mat"];
if isfile(spike_file)
    fprintf("Using existing binned spiking data from output/neuro/SpikeVar_spike_data.mat.\n")
    load(spike_file, 'spike_data_table');
else 
    fprintf("Extracting binned spiking data.\n")
    spike_data_table = a_prepare_spike_data(base_dir, data_dir, stim_dir);   
end

% 2. Calculate spike permutation entropy
%spike_data_table = b_compute_spike_pe(spike_data_table, neuro_dir); 

% alternatively load precomputed PE table
load([neuro_dir + 'SpikeVar_spike_data_PE.mat'], 'spike_data_table')

%% Stimuli: Compute feature estimates using HMAX and VGG
% load behavioral information and stimuli information
load([data_dir + 'behavior/SpikeVar_learn_recog_info.mat'], 'fulltable_learn', 'fulltable_recog');
stim_info = create_stim_info(image_dir, fullfile(base_dir, "output/")); 

% 1. Get HMAX estimates
hmax_file = [hmax_dir + 'HMAX_estimates.mat']; 
if isfile(hmax_file)
    fprintf("Using precomputed HMAX activation estimates.\n")
    load(hmax_file, 'HMAX_estimates');
else 
    fprintf("Estimating HMAX feature activations.\n")
    fprinft("This may take a couple of minutes.\n")
    HMAX_estimates = a_get_hmax_estimates(image_dir, hmax_dir);
end 

% 2. Run pca HMAX features estimates 
fprintf("Running PCA on HMAX feature estimates.\n")
hmax_pc1_table = b_hmax_pca(HMAX_estimates, stim_info, fulltable_learn, fulltable_recog, hmax_dir);

% 3. Get VGG estimates
% We're using the precomputed VGG16 estimates. To compute them again, run
% stimuli/c_get_vgg_estimates.py
fprintf("Running PCA on VGG16 feature estimates.\n")
vgg_pc1_table = d_vgg16_pca(vgg_dir, stim_info, fulltable_learn, fulltable_recog);


%% Partial Least Squares 
% Perform PLS analysis on hippocampus and amygdala data
hpc_pls_results = pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, 'hippocampus');
% amyg_pls_results = pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, 'amygdala');
