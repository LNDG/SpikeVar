clear
clc
base_dir = '/Users/kamp/PhD/spikevar/repo/SpikeVar';
stim_dir = fullfile(base_dir, 'stimuli/');
vgg_dir = fullfile(base_dir, 'output', 'vgg16/');
hmax_dir = fullfile(base_dir, 'output', 'hmax/');
data_dir = fullfile(base_dir, 'data/');
cd(stim_dir)

%  load info about recognition phase and new HMAX estimates
load([data_dir 'SpikeVar_learn_recog_info.mat'])
hmax_dir = '/Users/kamp/PhD/spikevar/repo/SpikeVar/output/hmax/';
load([hmax_dir  'HMAX_estims_allstims.mat'])

%% use only stimuli that have been used during encoding or recognition
stim_used = get_stim_used(info, fulltable_learn, fulltable_recog);
info(stim_used == 0,:) = [];

%% PCA on patch wise median
patch_wise_med_c1(stim_used==0,:) = [];
patch_wise_med_c2(stim_used==0,:) = [];
pca_med_c1 = get_layer_wise_pca(patch_wise_med_c1);
pca_med_c2 = get_layer_wise_pca(patch_wise_med_c2);

%% PCA on patch wise sum
patch_wise_sum_c1(stim_used==0,:) = [];
patch_wise_sum_c2(stim_used==0,:) = [];
pca_sum_c1 = get_layer_wise_pca(patch_wise_sum_c1);
pca_sum_c2 = get_layer_wise_pca(patch_wise_sum_c2);

%% same on patch wise spatial SD  
patch_wise_sd_c1(stim_used==0,:) = [];
patch_wise_sd_c2(stim_used==0,:) = [];
pca_sd_c1 = get_layer_wise_pca(patch_wise_sd_c1);
pca_sd_c2 = get_layer_wise_pca(patch_wise_sd_c2);


%% get across-patch SD to save later
% patch_sd_c1 = std(patch_wise_med_c1, [],2); 
% patch_sd_c2 = std(patch_wise_med_c2,[],2);

%% new table
% only have some basic info in there
pc1_table = table; pc1_table.participant = fulltable_learn.Participant;

%% get stimulus indices
stim_learn_idx = get_stim_learn_idx(info, fulltable_learn);

%% add new HMAX estimated to table
pc1_table.c1_med_score = pca_med_c1.score(stim_learn_idx,1);
pc1_table.c2_med_score = pca_med_c2.score(stim_learn_idx,1);

pc1_table.c1_sd_score = pca_sd_c1.score(stim_learn_idx,1);
pc1_table.c2_sd_score = pca_sd_c2.score(stim_learn_idx,1);

pc1_table.c1_sum_score = pca_sum_c1.score(stim_learn_idx,1);
pc1_table.c2_sum_score = pca_sum_c2.score(stim_learn_idx,1);

%pc1_table.c1_patch_sd = patch_sd_c1(stim_learn_idx,1);
%pc1_table.c2_patch_sd = patch_sd_c2(stim_learn_idx,1);

%% save estimates
save([hmax_dir 'HMAX_pca.mat'], 'pc1_table')
