clear
clc

% check if we are in the correct directory, change it if needed
base_dir = '/Users/kamp/PhD/spikevar/repo/SpikeVar';
stim_dir = fullfile(base_dir, 'stimuli/');
pls_dir = fullfile(base_dir, 'pls/');
hmax_dir = '/Users/kamp/PhD/spikevar/output/hmax_output/';
vgg_dir = '/Users/kamp/PhD/spikevar/output/vgg16_output/';
neuro_dir = '/Users/kamp/PhD/spikevar/output/neuro/';
cd(stim_dir)

% add toolboxes
addpath(genpath(fullfile(base_dir, 'toolboxes', 'pls_rank')));
% VGG PCs
load([vgg_dir, 'VGG16_pca.mat'], 'pc1_table');
vgg_pc1_table = pc1_table;
% HMAX PCs
load([hmax_dir 'HMAX_pca.mat'], 'pc1_table');
hmax_pc1_table = pc1_table;
% VGG PCs
load([vgg_dir, 'VGG16_all_encoding_recog_pc1_sd.mat'], 'pc1_table')
% spiking PE
load([neuro_dir 'SpikeVar_spike_data_PE.mat'], 'spike_data_table');
% get IDs
pat_ids = unique(spike_data_table.Participant);

%% loop across subjects
for pat_id = pat_ids'
    % get relevant indices
    % set of neurons (all, hippocampus <3  or amygdala >3)
    % 1 = right Hippocampius, 2 = left Hippocampus
    % 3 = right Amygdala, 4 = left Amygdala
    for i_ses = 1:2
        % select neurons in hippocampus (brainArea < 3)
        hpc_ids = find(spike_data_table.Participant==pat_id &...
                       spike_data_table.brainArea < 3 & ... 
                       spike_data_table.Session==i_ses);
        % only continue if we have data
        if isempty(hpc_ids)
            continue
        end
         
        % get number of trials (varies between participants and sessions)
        n_trials_all = unique(spike_data_table.NumTrials(spike_data_table.Participant==pat_id &...
             spike_data_table.brainArea<3 &spike_data_table.Session==i_ses));
        % save neuron info
        neuron_ids = unique(spike_data_table.NeuronCount(hpc_ids));
        all_neuron_ids = unique(spike_data_table.NeuronCount(hpc_ids));
        for i_neuron = 1:length(all_neuron_ids)
            all_neuron_area(i_neuron) = unique(spike_data_table.brainArea(...
                spike_data_table.NeuronCount == all_neuron_ids(i_neuron)));
        end
        % extract PE for different motif lengths
        pe_dat_2 = spike_data_table.spike_permEN_2(hpc_ids);
        pe_dat_3 = spike_data_table.spike_permEN_3(hpc_ids);
        pe_dat_4 = spike_data_table.spike_permEN_4(hpc_ids);
        
        % re-arrange data PE data
        n_neuron_all = length(neuron_ids);
        pe_dat_2 = reshape(pe_dat_2, n_trials_all, n_neuron_all);
        pe_dat_3 = reshape(pe_dat_3, n_trials_all, n_neuron_all);
        pe_dat_4 = reshape(pe_dat_4, n_trials_all, n_neuron_all);

        % drop neurons that are very silent
        drop_neurons = (mean(pe_dat_2)<.0001) | ...
                       (mean(pe_dat_3)<.0001) | ... 
                       (mean(pe_dat_4)<.0001);
        pe_dat_2(:, drop_neurons) = [];
        pe_dat_3(:, drop_neurons) = [];
        pe_dat_4(:, drop_neurons) = [];
        used_neuron_ids = all_neurons_ids;
        used_neuron_ids(drop_neurons) = [];
        n_neuron_used = length(all_neuron_ids);
        
        % also save average PE of remaining neurons
        avg_pe_mat_2 = mean(pe_dat_2);
        avg_pe_mat_3 = mean(pe_dat_3);
        avg_pe_mat_4 = mean(pe_dat_4);

        % another check if we still have neurons left
        if size(pe_dat_2,2) == 0
            continue
        end
               
        % get trials in which less than 1/3 of neurons are
        % spiking at all
        drop_trials = find(sum(sign(pe_dat_2),2)<=...
            (size(pe_dat_2,2)/3));        
        
        % initialize the spiking input for the pls and drop silent trials
        neuro_data = [pe_dat_2, pe_dat_3, pe_dat_4];
        neuro_data(drop_trials,:) = [];
        n_trails_used = size(neuro_data,1);
        
        % initialize stimulus input for pls
        stim_data = struct();
        
        % prepare VGG data
        layers = {'first', 'second', 'third', 'fourth', 'fifth'};
        for layer = layers
            % vgg sd 
            col_name = [layer '_layer_sd_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(vgg_pc1_table(:, col_name), ... 
                hpc_ids, n_trials_all, n_neuron_all, drop_trials);        
            % vgg sum 
            col_name = [layer '_layer_sum_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(vgg_pc1_table(:, col_name), ... 
                hpc_ids, n_trials_all, n_neuron_all, drop_trials); 
        end
        % prepare hmax data
        layers = {'c1', 'c2'};
        for layer = layers
            % hmax sd 
            col_name = [layer '_sd_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(hmax_pc1_table(:, col_name), ... 
                hpc_ids, n_trials_all, n_neuron_all, drop_trials); 
            % hmax sd 
            col_name = [layer '_sum_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(hmax_pc1_table(:, col_name), ... 
                hpc_ids, n_trials_all, n_neuron_all, drop_trials);
        end
           
        % run PLS models
        % set up some more options for PLS
        option = [];
        option.method = 3; % behavioural PLS
        option.num_perm = 1000; % number of permutations
        option.num_boot = 1000; % number of bootstraps
        option.cormode = 8; % rank correlation

        % set up the stimulus info matrix for 
        % a model that includes all vgg and HMAX layers
        option.stacked_behavdata = [
            stim_data.first_layer_sum,...                                  % vgg sum
            stim_data.second_layer_sum,...
            stim_data.third_layer_sum,...
            stim_data.fourth_layer_sum,...
            stim_data.fifth_layer_sum,...
            stim_data.first_layer_sd,...                                   % vgg sd
            stim_data.second_layer_sd,...
            stim_data.third_layer_sd,...
            stim_data.fourth_layer_sd,...
            stim_data.fifth_layer_sd,...
            stim_data.c1_sum,...                                           % hmax sum
            stim_data.c2_sum,...
            stim_data.c1_sd,...                                            % hmax sd
            stim_data.c2_sd];

        % run model
        % neural matrix consists of neuron and trial-wise PE for motif lengths 2, 3, & 4
        result_combined = pls_analysis({neuro_data}, n_trails_used, 1, option);
       
        % run one model using late layers
        option.stacked_behavdata = [
            stim_data.third_layer_sum,...
            stim_data.fourth_layer_sum,...
            stim_data.fifth_layer_sum,...
            stim_data.third_layer_sd,...
            stim_data.fourth_layer_sd,...
            stim_data.fifth_layer_sd,...
            stim_data.c2_sum,...
            stim_data.c2_sd];
        result_late = pls_analysis({neuro_data}, n_trails_used, 1, option);

        % run one model using early layers
        option.stacked_behavdata = [
            stim_data.first_layer_sum,...                                  
            stim_data.second_layer_sum,...
            stim_data.third_layer_sum,...
            stim_data.first_layer_sd,...                                   
            stim_data.second_layer_sd,...
            stim_data.third_layer_sd,...
            stim_data.c1_sum,...                                           
            stim_data.c1_sd];                                               
        result_early = pls_analysis(datamat_lst, n_trails_used, 1, option);

        % save subject-specific stimulus weights
        behav_weights.all_layers{pat_id,i_ses} = result_combined.v;
        behav_weights.late_layers{pat_id,i_ses} = result_late.v;
        behav_weights.early_layers{pat_id,i_ses} = result_early.v;

        % save number of trials
        trail_nums(pat_id, i_ses) = n_trails_used;
        % number of neurons
        neuron_nums(pat_id, i_ses) = n_neurons_used;

        % save LVLV corrs
        lvlv_corrs.all_layers{i_ses}(pat_id) = result_combined.lvlvcorr(1);
        lvlv_corrs.late_layers{i_ses}(pat_id) = result_late.lvlvcorr(1);
        lvlv_corrs.early_layers{i_ses}(pat_id) = result_early.lvlvcorr(1);

        % also get bounds of LVLV corrs  
        % use the actual lvlv corr distribution
        lvlv_corr_lims.all_layers{i_ses}(pat_id,:) = [percentile(squeeze(...
            result_combined.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_combined.boot_result.lvlvcorrdistrib(1,:)),99)];
        lvlv_corr_lims.late_layers{i_ses}(pat_id,:) = [percentile(squeeze(...
            result_late.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_late.boot_result.lvlvcorrdistrib(1,:)),99)];
        lvlv_corr_lims.early_layers{i_ses}(pat_id,:) = [percentile(squeeze(...
            result_early.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_early.boot_result.lvlvcorrdistrib(1,:)),99)];
    end
end
% save stuff
save('SpikeVar_early_vs_late_PE_PLS_results.mat', 'behav_weights',...
    'lvlv_corrs','lvlv_corr_lims', 'trial_nums', 'neuron_nums')
