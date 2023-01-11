clear
clc

% check if we are in the correct directory, change it if needed
cd('/Users/waschke/OneDrive/Old data and code/SV/OurCode/')
%  load info about recognition phase
load('SpikeVar_learn_recog_info.mat')
% also load info about recognition phase
load('SPSSinput_learn_recog.mat')
% also load new HMAX estimates
hmax_dir = '/Users/waschke/OneDrive/Old data and code/HMAX output/';
vggdir = '/Users/waschke/OneDrive/Old data and code/VGG16_output/';
load([hmax_dir  'HMAX_estims_allstims_April21_22.mat'])
% load task variant information
task_vars = import_task_variants([vggdir 'SpikeVar_stim_task_variants.csv']);
%% which stimuli were presented during the encoding phase?
for s_i = 1:size(all_HMAX_info,1)
   stim_use_learn(s_i) = ismember(all_HMAX_info.StimulusCode(s_i), fulltable_learn.stimuliLearn) ;
   stim_use_recog(s_i) = ismember(all_HMAX_info.StimulusCode(s_i), fulltable_recog.stimuliRecog);
end
all_stim_use  = stim_use_learn+stim_use_recog;
%% also load VGG16 activation estimates
VGG16lay1feats = importfile_VGG_feat1([vggdir 'SpikeVar_VGG16_lay1_feature_sum.csv']);
% get rid of first column
VGG16lay1feats(:,1) = [];
% make it an array
VGG16lay1feats = table2array(VGG16lay1feats);
% replace NaNs with zeros
VGG16lay1feats(isnan(VGG16lay1feats)) = 0;
% kick out non used stimuli
VGG16lay1feats(all_stim_use==0,:) = [];

VGG16lay2feats = importfile_VGG_feat2([vggdir 'SpikeVar_VGG16_lay2_feature_sum.csv']);
% get rid of first column
VGG16lay2feats(:,1) = [];
% make it an array
VGG16lay2feats = table2array(VGG16lay2feats);
% replace NaNs with zeros
VGG16lay2feats(isnan(VGG16lay2feats)) = 0;
% kick out non used stimuli
VGG16lay2feats(all_stim_use==0,:) = [];

VGG16lay3feats = importfile_VGG_feat3([vggdir 'SpikeVar_VGG16_lay3_feature_sum.csv']);
% get rid of first column
VGG16lay3feats(:,1) = [];
% make it an array
VGG16lay3feats = table2array(VGG16lay3feats);
% replace NaNs with zeros
VGG16lay3feats(isnan(VGG16lay3feats)) = 0;
VGG16lay3feats(all_stim_use==0,:) = [];

VGG16lay4feats = importfile_VGG_feat4([vggdir 'SpikeVar_VGG16_lay4_feature_sum.csv']);
% get rid of first column
VGG16lay4feats(:,1) = [];
% make it an array
VGG16lay4feats = table2array(VGG16lay4feats);
% replace NaNs with zeros
VGG16lay4feats(isnan(VGG16lay4feats)) = 0;
VGG16lay4feats(all_stim_use==0,:) = [];

VGG16lay5feats = importfile_VGG_feat5([vggdir 'SpikeVar_VGG16_lay5_feature_sum.csv']);
% get rid of first column
VGG16lay5feats(:,1) = [];
% make it an array
VGG16lay5feats = table2array(VGG16lay5feats);
% replace NaNs with zeros
VGG16lay5feats(isnan(VGG16lay5feats)) = 0;
VGG16lay5feats(all_stim_use==0,:) = [];


VGG16layffeats = importfile_VGG_feat_final([vggdir 'SpikeVar_VGG16_lay_final_features.csv']);
% get rid of first column
VGG16layffeats(:,1) = [];
VGG16layffeats(1,:) = [];
% make it an array
VGG16layffeats = table2array(VGG16layffeats);
% replace NaNs with zeros
VGG16layffeats(isnan(VGG16layffeats)) = 0;
VGG16layffeats(isinf(VGG16layffeats)) = 0;
% kick out completely empty features
VGG16layffeats(all_stim_use==0,:) = [];
% kick out completely empty features
VGG16layffeats_sum = sum(VGG16layffeats,1);
VGG16layffeats(:,VGG16layffeats_sum==0) = [];

% also clean task_variant vector
task_vars(all_stim_use==0,:) = [];
all_HMAX_info(all_stim_use ==0,:) = [];

%% run PCA on ech feature map

[coeff.vgg1,score.vgg1,...
    latent.vgg1,~,exp_var.vgg1,~] =...
    pca(zscore(VGG16lay1feats(:, [1:4,6:end])), 'Centered', true,'VariableWeights', 'Variance');
loadings.vgg1 = (corr(score.vgg1,...
    VGG16lay1feats(:,1:size(score.vgg1,2))))';
exp_var_1st_pc.vgg1 = exp_var.vgg1(1);
exp_normvar_1st_pc.vgg1 = exp_var.vgg1(1)/...
                    (100/size(exp_var.vgg1,1));

[coeff.vgg2,score.vgg2,...
    latent.vgg2,~,exp_var.vgg2,~] =...
    pca(zscore(VGG16lay2feats), 'Centered', true,'VariableWeights', 'Variance');
loadings.vgg2 = (corr(score.vgg2,...
    VGG16lay2feats(:,1:size(score.vgg2,2))))';
exp_var_1st_pc.vgg2 = exp_var.vgg2(1);
exp_normvar_1st_pc.vgg2 = exp_var.vgg2(1)/...
                    (100/size(exp_var.vgg2,1));


[coeff.vgg3,score.vgg3,...
    latent.vgg3,~,exp_var.vgg3,~] =...
    pca(zscore(VGG16lay3feats), 'Centered', true,'VariableWeights', 'Variance');
loadings.vgg3 = (corr(score.vgg3,...
    VGG16lay3feats(:,1:size(score.vgg3,2))))';
exp_var_1st_pc.vgg3 = exp_var.vgg3(1);
exp_normvar_1st_pc.vgg3 = exp_var.vgg3(1)/...
                    (100/size(exp_var.vgg3,1));
                
[coeff.vgg4,score.vgg4,...
    latent.vgg4,~,exp_var.vgg4,~] =...
    pca(zscore(VGG16lay4feats), 'Centered', true,'VariableWeights', 'Variance');
loadings.vgg4 = (corr(score.vgg4,...
    VGG16lay4feats(:,1:size(score.vgg4,2))))';
exp_var_1st_pc.vgg4 = exp_var.vgg4(1);
exp_normvar_1st_pc.vgg4 = exp_var.vgg4(1)/...
                    (100/size(exp_var.vgg4,1));

[coeff.vgg5,score.vgg5,...
    latent.vgg5,~,exp_var.vgg5,~] =...
    pca(zscore(VGG16lay5feats), 'Centered', true,'VariableWeights', 'Variance');
loadings.vgg5 = (corr(score.vgg5,...
    VGG16lay5feats(:,1:size(score.vgg5,2))))';
exp_var_1st_pc.vgg5 = exp_var.vgg5(1);
exp_normvar_1st_pc.vgg5 = exp_var.vgg5(1)/...
                    (100/size(exp_var.vgg5,1));
                
                
[coeff.vggf,score.vggf,...
    latent.vggf,~,exp_var.vggf,~] =...
    pca(zscore(VGG16layffeats), 'Centered', true,'VariableWeights', 'Variance');
loadings.vggf = (corr(score.vggf,...
    VGG16layffeats(:,1:size(score.vggf,2))))';
exp_var_1st_pc.vggf = exp_var.vggf(1);
exp_normvar_1st_pc.vggf = exp_var.vggf(1)/...
                    (100/size(exp_var.vggf,1));
                
%% determine highest variance explained to get the nmber of components needed for each layer
max_1st_pc_var = max(struct2array(exp_var_1st_pc));

num_comp_need.vgg1 = find(cumsum(exp_var.vgg1)>=max_1st_pc_var,1); 
num_comp_need.vgg2 = find(cumsum(exp_var.vgg2)>=max_1st_pc_var,1); 
num_comp_need.vgg3 = find(cumsum(exp_var.vgg3)>=max_1st_pc_var,1); 
num_comp_need.vgg4 = find(cumsum(exp_var.vgg4)>=max_1st_pc_var,1); 
num_comp_need.vgg5 = find(cumsum(exp_var.vgg5)>=max_1st_pc_var,1); 
%% new table
% only have some basic info in there
VGG1_sum_pc1_table = table; VGG1_sum_pc1_table.Participant = fulltable.Participant;
VGG2_sum_pc1_table = table; VGG2_sum_pc1_table.Participant = fulltable.Participant;
VGG3_sum_pc1_table = table; VGG3_sum_pc1_table.Participant = fulltable.Participant;
VGG4_sum_pc1_table = table; VGG4_sum_pc1_table.Participant = fulltable.Participant;
VGG5_sum_pc1_table = table; VGG5_sum_pc1_table.Participant = fulltable.Participant;
VGGf_sum_pc1_table = table; VGGf_sum_pc1_table.Participant = fulltable.Participant;
%% average scores of pcs needed to explain 30 % of variance
VGG1_sum_pcavg_table = table; VGG1_sum_pcavg_table.Participant = fulltable.Participant;
VGG2_sum_pcavg_table = table; VGG2_sum_pcavg_table.Participant = fulltable.Participant;
VGG3_sum_pcavg_table = table; VGG3_sum_pcavg_table.Participant = fulltable.Participant;
VGG4_sum_pcavg_table = table; VGG4_sum_pcavg_table.Participant = fulltable.Participant;
VGG5_sum_pcavg_table = table; VGG5_sum_pcavg_table.Participant = fulltable.Participant;
%% instead use the first 3 components
VGG1_sum_pc_table = table; VGG1_sum_pc_table.Participant = fulltable.Participant;
VGG2_sum_pc_table = table; VGG2_sum_pc_table.Participant = fulltable.Participant;
VGG3_sum_pc_table = table; VGG3_sum_pc_table.Participant = fulltable.Participant;
VGG4_sum_pc_table = table; VGG4_sum_pc_table.Participant = fulltable.Participant;
VGG5_sum_pc_table = table; VGG5_sum_pc_table.Participant = fulltable.Participant;
%% get stimulus indices
for il = 1:size(fulltable,1)
    stimind_f(il) = find(all_HMAX_info.StimulusCode==fulltable.stimuliLearn(il));
end
%% add new HMAX estimated to table
VGG1_sum_pc1_table.lay1_score = score.vgg1(stimind_f,1);
VGG2_sum_pc1_table.lay2_score = score.vgg2(stimind_f,1);
VGG3_sum_pc1_table.lay3_score = score.vgg3(stimind_f,1);
VGG4_sum_pc1_table.lay4_score = score.vgg4(stimind_f,1);
VGG5_sum_pc1_table.lay5_score = score.vgg5(stimind_f,1);
VGGf_sum_pc1_table.layf_score = score.vggf(stimind_f,1);

VGG1_sum_pcavg_table.lay1_score = score.vgg1(stimind_f,1);
VGG2_sum_pcavg_table.lay2_score = score.vgg2(stimind_f,1);
VGG3_sum_pcavg_table.lay3_score = mean(score.vgg3(stimind_f,1:2),2);
VGG4_sum_pcavg_table.lay4_score = mean(score.vgg4(stimind_f,1:3),2);
VGG5_sum_pcavg_table.lay5_score = mean(score.vgg5(stimind_f,1:6),2);

VGG1_sum_pc_table.lay1_score = score.vgg1(stimind_f,1:num_comp_need.vgg1);
VGG2_sum_pc_table.lay2_score = score.vgg2(stimind_f,1:num_comp_need.vgg2);
VGG3_sum_pc_table.lay3_score = score.vgg3(stimind_f,1:num_comp_need.vgg3);
VGG4_sum_pc_table.lay4_score = score.vgg4(stimind_f,1:num_comp_need.vgg4);
VGG5_sum_pc_table.lay5_score = score.vgg5(stimind_f,1:num_comp_need.vgg5);


%% also add features
VGG1_sum_pc1_table.sum_feats = VGG16lay1feats(stimind_f,:);
VGG2_sum_pc1_table.sum_feats = VGG16lay2feats(stimind_f,:);
VGG3_sum_pc1_table.sum_feats = VGG16lay3feats(stimind_f,:);
VGG4_sum_pc1_table.sum_feats = VGG16lay4feats(stimind_f,:);
VGG5_sum_pc1_table.sum_feats = VGG16lay5feats(stimind_f,:);
VGGf_sum_pc1_table.sum_feats = VGG16layffeats(stimind_f,:);
%% save estimates
save([vggdir 'VGG16_all_encoding_recog_layer1_sum_pc1.mat'], 'VGG1_sum_pc1_table', '-v7.3')
save([vggdir 'VGG16_all_encoding_recog_layer2_sum_pc1.mat'], 'VGG2_sum_pc1_table', '-v7.3')
save([vggdir 'VGG16_all_encoding_recog_layer3_sum_pc1.mat'], 'VGG3_sum_pc1_table', '-v7.3')
save([vggdir 'VGG16_all_encoding_recog_layer4_sum_pc1.mat'], 'VGG4_sum_pc1_table', '-v7.3')
save([vggdir 'VGG16_all_encoding_recog_layer5_sum_pc1.mat'], 'VGG5_sum_pc1_table', '-v7.3')
save([vggdir 'VGG16_all_encoding_recog_layerf_sum_pc1.mat'], 'VGGf_sum_pc1_table', '-v7.3')

save([vggdir 'VGG16_all_encoding_recog_layer1_sum_pcavg.mat'], 'VGG1_sum_pcavg_table')
save([vggdir 'VGG16_all_encoding_recog_layer2_sum_pcavg.mat'], 'VGG2_sum_pcavg_table')
save([vggdir 'VGG16_all_encoding_recog_layer3_sum_pcavg.mat'], 'VGG3_sum_pcavg_table')
save([vggdir 'VGG16_all_encoding_recog_layer4_sum_pcavg.mat'], 'VGG4_sum_pcavg_table')
save([vggdir 'VGG16_all_encoding_recog_layer5_sum_pcavg.mat'], 'VGG5_sum_pcavg_table')

save([vggdir 'VGG16_all_encoding_recog_layer1_sum_pcs.mat'], 'VGG1_sum_pc_table')
save([vggdir 'VGG16_all_encoding_recog_layer2_sum_pcs.mat'], 'VGG2_sum_pc_table')
save([vggdir 'VGG16_all_encoding_recog_layer3_sum_pcs.mat'], 'VGG3_sum_pc_table')
save([vggdir 'VGG16_all_encoding_recog_layer4_sum_pcs.mat'], 'VGG4_sum_pc_table')
save([vggdir 'VGG16_all_encoding_recog_layer5_sum_pcs.mat'], 'VGG5_sum_pc_table')

%% correlate scores across layers
comb_table = table(score.vgg1(:,1),score.vgg2(:,1),score.vgg3(:,1), score.vgg4(:,1),...
    score.vgg5(:,1),...
    'Variablenames', {'VGG1','VGG2', 'VGG3', 'VGG4', 'VGG5'});

figure
corrplot(comb_table,'Type',"Pearson",'TestR',"on")
