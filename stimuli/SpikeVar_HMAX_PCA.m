clear
clc
% check if we are in the correct directory, change it if needed
cd('/Users/kamp/PhD/spikevar/repo/SpikeVar/stimuli/')
%  load info about recognition phase
load('SpikeVar_learn_recog_info.mat')
% also load new HMAX estimates
hmax_dir = '/Users/kamp/PhD/spikevar/hmax_output/';
load([hmax_dir  'HMAX_estims_allstims.mat'])

%% which stimuli were presented during the encoding and recognition phase?
for s_i = 1:size(all_HMAX_info,1)
   stim_use_learn(s_i) = ismember(all_HMAX_info.StimulusCode(s_i), fulltable_learn.stimuliLearn) ;
   stim_use_recog(s_i) = ismember(all_HMAX_info.StimulusCode(s_i), fulltable_recog.stimuliRecog);
end
all_stim_use  = stim_use_learn+stim_use_recog;
all_HMAX_info(all_stim_use ==0,:) = [];
% clean stimulus info vectors
patch_wise_med_c1(all_stim_use ==0,:) = []; % why not for patch_wise_sd_c1
patch_wise_med_c2(all_stim_use ==0,:) = []; 

%% PCA on patch-wise C1 an C2 estimates
pca_med_c1 = get_patch_wise_pca(patch_wise_med_c1);
pca_med_c2 = get_patch_wise_pca(patch_wise_med_c2);
                
%% same on patch wise spatial SD               
pca_sd_c1 = get_patch_wise_pca(patch_wise_sd_c1);
pca_sd_c2 = get_patch_wise_pca(patch_wise_sd_c2);
  
                                               
%% get across-patch SD to save later
c1_patch_sd = std(patch_wise_med_c1, [],2); 
c2_patch_sd = std(patch_wise_med_c2,[],2);

%% new table
% only have some basic info in there
C1_pc1_table = table; C1_pc1_table.Participant = fulltable.Participant;
C2_pc1_table = table; C2_pc1_table.Participant = fulltable.Participant;

%% get stimulus indices
for il = 1:size(fulltable,1)
    stimind_f(il) = find(all_HMAX_info.StimulusCode==fulltable.stimuliLearn(il)); % stimind_f not defined
end
%% add new HMAX estimated to table
C1_pc1_table.C1_med_score = pca_med_c1.score(stimind_f,1);
C2_pc1_table.C2_med_score = pca_med_c2.score(stimind_f,1);

C1_pc1_table.C1_sd_score = pca_sd_c1.score(stimind_f,1);
C2_pc1_table.C2_sd_score = pca_sd_c2.score(stimind_f,1);

C1_pc1_table.C1_patch_sd = c1_patch_sd(stimind_f,1);
C2_pc1_table.C2_patch_sd = c2_patch_sd(stimind_f,1);

%% save estimates
save([hmax_dir 'HMAX_all_encoding_recog_C1_pc1.mat'], 'C1_pc1_table')
save([hmax_dir 'HMAX_all_encoding_recog_C2_pc1.mat'], 'C2_pc1_table')


