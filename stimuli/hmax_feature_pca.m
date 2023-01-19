clear
clc
% check if we are in the correct directory, change it if needed
base_dir = '/Users/kamp/PhD/spikevar/repo/SpikeVar/';
stim_dir = [base_dir 'stimuli/'];
cd(stim_dir)
%  load info about recognition phase and new HMAX estimates
load('SpikeVar_learn_recog_info.mat')
hmax_dir = '/Users/kamp/PhD/spikevar/output/hmax_output/';
load([hmax_dir  'HMAX_estims_allstims.mat'])

%% PCA on patch-wise C1 an C2 estimates
pca_med_c1 = get_layer_wise_pca(patch_wise_med_c1);
pca_med_c2 = get_layer_wise_pca(patch_wise_med_c2);
                
%% same on patch wise spatial SD               
pca_sd_c1 = get_layer_wise_pca(patch_wise_sd_c1);
pca_sd_c2 = get_layer_wise_pca(patch_wise_sd_c2);
                                                 
%% get across-patch SD to save later
patch_sd_c1 = std(patch_wise_med_c1, [],2); 
patch_sd_c2 = std(patch_wise_med_c2,[],2);

%% new table
% only have some basic info in there
pc1_table = table; c1_pc1_table.participant = fulltable_learn.Participant;

%% get stimulus indices
stim_learn_idx = get_stim_learn_idx(all_HMAX_info, fulltable_learn);

%% add new HMAX estimated to table
pc1_table.c1_med_score = pca_med_c1.score(stim_learn_idx,1);
pc1_table.c2_med_score = pca_med_c2.score(stim_learn_idx,1);

pc1_table.c1_sd_score = pca_sd_c1.score(stim_learn_idx,1);
pc1_table.c2_sd_score = pca_sd_c2.score(stim_learn_idx,1);

pc1_table.c1_patch_sd = patch_sd_c1(stim_learn_idx,1);
pc1_table.c2_patch_sd = patch_sd_c2(stim_learn_idx,1);

%% save estimates
save([hmax_dir 'HMAX_pca.mat'], 'pc1_table')
