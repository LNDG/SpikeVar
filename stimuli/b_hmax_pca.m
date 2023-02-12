function pc1_table = b_hmax_pca(HMAX_estimates, stim_info, fulltable_learn, fulltable_recog, hmax_dir)
%b_hmax_pca computes the pca across the extraceted HMAX features
%   Takes as input the HMAX estimates, the info table, and the
%   tables containing trialwise stimulus information

patch_wise_med_c1 = HMAX_estimates.patch_wise_med_c1;
patch_wise_med_c2 = HMAX_estimates.patch_wise_med_c2;
patch_wise_sum_c1 = HMAX_estimates.patch_wise_sum_c1;
patch_wise_sum_c2 = HMAX_estimates.patch_wise_sum_c2;
patch_wise_sd_c1 = HMAX_estimates.patch_wise_sd_c1;
patch_wise_sd_c2 = HMAX_estimates.patch_wise_sd_c2;

%% use only stimuli that have been used during encoding or recognition
stim_used = get_stim_used(stim_info, fulltable_learn, fulltable_recog);
stim_info(stim_used == 0,:) = [];

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

%% new table
% only have some basic info in there
pc1_table = table; pc1_table.participant = fulltable_learn.Participant;

%% get stimulus indices
stim_learn_idx = get_stim_learn_idx(stim_info, fulltable_learn);

%% add new HMAX estimated to table
pc1_table.c1_med_score = pca_med_c1.score(stim_learn_idx,1);
pc1_table.c2_med_score = pca_med_c2.score(stim_learn_idx,1);

pc1_table.c1_sd_score = pca_sd_c1.score(stim_learn_idx,1);
pc1_table.c2_sd_score = pca_sd_c2.score(stim_learn_idx,1);

pc1_table.c1_sum_score = pca_sum_c1.score(stim_learn_idx,1);
pc1_table.c2_sum_score = pca_sum_c2.score(stim_learn_idx,1);

%% save estimates
save([hmax_dir + 'HMAX_pca.mat'], 'pc1_table')

end

