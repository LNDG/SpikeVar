clear
clc
base_dir = '/Users/kamp/PhD/spikevar/repo/SpikeVar';
stim_dir = fullfile(base_dir, 'stimuli/');
vgg_dir = fullfile(base_dir, 'output', 'vgg16/');
hmax_dir = fullfile(base_dir, 'output', 'hmax/');
data_dir = fullfile(base_dir, 'data/');
cd(stim_dir)

%  load info about recognition phase and new HMAX estimates
load([data_dir 'SpikeVar_learn_recog_info.mat']);
load([hmax_dir 'HMAX_estims_allstims.mat']);

% load task variant information
vgg_info = readtable([vgg_dir 'SpikeVar_VGG16_Complexity.csv']);
% layer_names 
layer_names = {'first', 'second', 'third', 'fourth', 'fifth'};

%% VGG16 feature
% get stimuli that have been used during encoding/recognition
stim_used = get_stim_used(info, fulltable_learn, fulltable_recog);
info(stim_used ==0,:) = [];
% get the stimulus indexes in the learned order for each participant
stim_learn_idx = get_stim_learn_idx(info, fulltable_learn);

% initiate output tables
pc1_table = table; pc1_table.participant = fulltable_learn.Participant;
pc_table = table; pc_table.participant = fulltable_learn.Participant;
% table to store averages scores across PCs to get to 30% explained variance
pc_avg_table = table; pc_avg_table.participant = fulltable_learn.Participant;

stats = {'sum', 'nz_sd', 'nz_num'};
for i = 1:length(stats)
    max_1st_pc_var = 0;
    for j = 1:length(layer_names)
        layer = [layer_names{j}, '_layer'];
        % load layer-wise feature values
        file_path = [vgg_dir, 'SpikeVar_VGG16_', layer_names{j}, '_mp-layer_features_', stats{i}, '.csv'];
        layer_wise_features = load_vgg16_layer_features(file_path, stim_used);
        if j==1
            layer_wise_features(:,5) = []; % feature with zero variance
        end
        % apply pca
        layer_wise_pca.(layer) = get_layer_wise_pca(layer_wise_features);
        % determine highest variance explained explained to get num of components for
        % each layer
        max_1st_pc_var = max([max_1st_pc_var, layer_wise_pca.(layer).exp_var_1st_pc]);
    end
      
    % save results
    for j = 1:length(layer_names)
        layer = [layer_names{j}, '_layer'];
        num_comp_need = find(cumsum(layer_wise_pca.(layer).exp_var)>=max_1st_pc_var,1);
        % add scores to tables
        col_name = [layer_names{j} '_layer_' stats{i} '_score'];
        layer_pca_score = layer_wise_pca.(layer).score;
        pc1_table.(col_name) = layer_pca_score(stim_learn_idx,1);
        pca_table = layer_pca_score;
        pc_table.(col_name) = layer_pca_score(stim_learn_idx,1:num_comp_need);
        if j <= 3 
            pc_avg_table.(col_name) = layer_pca_score(stim_learn_idx,1);
        elseif j == 4
            pc_avg_table.(col_name) = mean(layer_pca_score(stim_learn_idx,1:3),2); % inconsistent between sd and sum pca
        elseif j == 5
            pc_avg_table.(col_name) = mean(layer_pca_score(stim_learn_idx,1:8),2);
        end
    end
end 
save([vgg_dir 'VGG16_pca.mat'], 'pc1_table', 'pc_table', 'pc_avg_table', 'pca_table')

