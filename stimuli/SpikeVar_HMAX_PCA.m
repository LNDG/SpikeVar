clear
clc
% check if we are in the correct directory, change it if needed
cd('/Users/waschke/OneDrive/Old data and code/SV/OurCode/')
%  load info about recognition phase
load('SpikeVar_learn_recog_info.mat')
% also load new HMAX estimates
hmax_dir = '/Users/waschke/OneDrive/Old data and code/HMAX output/';
load([hmax_dir  'HMAX_estims_allstims_April21_22.mat'])
vggdir = '/Users/waschke/OneDrive/Old data and code/VGG16_output/';
% load task variant information
task_vars = import_task_variants([vggdir 'SpikeVar_stim_task_variants.csv']);
%% which stimuli were presented during the encoding and recognition phase?
for s_i = 1:size(all_HMAX_info,1)
   stim_use_learn(s_i) = ismember(all_HMAX_info.StimulusCode(s_i), fulltable_learn.stimuliLearn) ;
   stim_use_recog(s_i) = ismember(all_HMAX_info.StimulusCode(s_i), fulltable_recog.stimuliRecog);
end
all_stim_use  = stim_use_learn+stim_use_recog;
%  clean task_variant vector
task_vars(all_stim_use==0,:) = [];
all_HMAX_info(all_stim_use ==0,:) = [];
% clean stimulus info vectors
patch_wise_c1(all_stim_use ==0,:) = [];
patch_wise_c2(all_stim_use ==0,:) = [];
%% PCA on patch-wise C1 an C2 estimates
[coeff.c1,score.c1,...
    latent.c1,~,exp_var.c1,~] =...
    pca(zscore(patch_wise_c1), 'Centered', true,'VariableWeights', 'Variance');

loadings.c1 = (corr(score.c1,...
    patch_wise_c1(:,1:size(score.c1,2))))';
exp_var_1st_pc.c1 = exp_var.c1(1);
exp_normvar_1st_pc.c1 = exp_var.c1(1)/...
                    (100/size(exp_var.c1,1));


[coeff.c2,score.c2,...
    latent.c2,~,exp_var.c2,~] =...
    pca(zscore(patch_wise_c2), 'Centered', true,'VariableWeights', 'Variance');
loadings.c2= (corr(score.c2,...
    patch_wise_c2(:,1:size(score.c2,2))))';
exp_var_1st_pc.c2 = exp_var.c2(1);
exp_normvar_1st_pc.c2 = exp_var.c2(1)/...
                    (100/size(exp_var.c2,1));
                
%% same on patch wise spatial SD               
[coeff.sd_c1,score.sd_c1,...
    latent.sd_c1,~,exp_var.sd_c1,~] =...
    pca(zscore(patch_wise_c1sd), 'Centered', true,'VariableWeights', 'Variance');
exp_var_1st_pc.sd_c1 = exp_var.sd_c1(1);
exp_normvar_1st_pc.sd_c1 = exp_var.sd_c1(1)/...
                    (100/size(exp_var.sd_c1,1));


[coeff.sd_c2,score.sd_c2,...
    latent.sd_c2,~,exp_var.sd_c2,~] =...
    pca(zscore(patch_wise_c2sd), 'Centered', true,'VariableWeights', 'Variance');
loadings.sd_c2= (corr(score.sd_c2,...
    patch_wise_c2sd(:,1:size(score.sd_c2,2))))';
exp_var_1st_pc.sd_c2 = exp_var.sd_c2(1);
exp_normvar_1st_pc.sd_c2 = exp_var.sd_c2(1)/...
                    (100/size(exp_var.sd_c2,1));    
                                               
%% also quickly get across-patch SD to save later
c1_patch_sd = std(patch_wise_c1, [],2);
c2_patch_sd = std(patch_wise_c2,[],2);
%% new table
% only have some basic info in there
C1_pc1_table = table; C1_pc1_table.Participant = fulltable.Participant;
C2_pc1_table = table; C2_pc1_table.Participant = fulltable.Participant;

%% get stimulus indices
for il = 1:size(fulltable,1)
    stimind_f(il) = find(all_HMAX_info.StimulusCode==fulltable.stimuliLearn(il));
end
%% add new HMAX estimated to table
C1_pc1_table.C1_score = score.c1(stimind_f,1);
C2_pc1_table.C2_score = score.c2(stimind_f,1);

C1_pc1_table.C1_sd_score = score.sd_c1(stimind_f,1);
C2_pc1_table.C2_sd_score = score.sd_c2(stimind_f,1);

C1_pc1_table.C1_patch_sd = c1_patch_sd(stimind_f,1);
C2_pc1_table.C2_patch_sd = c2_patch_sd(stimind_f,1);

%% save estimates
save([hmax_dir 'HMAX_all_encoding_recog_C1_pc1.mat'], 'C1_pc1_table')
save([hmax_dir 'HMAX_all_encoding_recog_C2_pc1.mat'], 'C2_pc1_table')


