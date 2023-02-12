function pls_results = pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, region)
%pls_subject_level Function that loops over all subjects and performs 3 behavioral pls
%   on for each subject.
%   Takes as input the spike_data_table, and the first pca of hmax and vgg 
%   outputs.

%% loop across subjects
count = 1; 
% get IDs
pat_ids = unique(spike_data_table.Participant);

for pat_id = pat_ids'
    for i_ses = 1:2
        % get relevant indices
        % set of neurons (all, hippocampus <3  or amygdala =>3)
        % 1 = right Hippocampius, 2 = left Hippocampus
        % 3 = right Amygdala, 4 = left Amygdala
        if region == "hippocampus"
            % select neurons in hippocampus (brainArea < 3)
            area_ids = find(spike_data_table.Participant==pat_id &...
                           spike_data_table.brainArea < 3 & ... 
                           spike_data_table.Session==i_ses);
        elseif region == "amygdala"
            % select neurons in amygdala (brainArea => 3)
            area_ids = find(spike_data_table.Participant==pat_id &...
                           spike_data_table.brainArea >= 3 & ... 
                           spike_data_table.Session==i_ses)
        else 
            fprintf('Region has to hippocampus or amygdala')
        end
        
        % only continue if we have data
        if isempty(area_ids)
            continue
        end
     
        % get number of trials (varies between participants and sessions)
        n_trials_all = unique(spike_data_table.NumTrials(spike_data_table.Participant==pat_id &...
             spike_data_table.brainArea<3 &spike_data_table.Session==i_ses));
        % save neuron info
        neuron_ids = unique(spike_data_table.NeuronCount(area_ids));

        % extract PE for different motif lengths
        pe_dat_2 = spike_data_table.spike_permEn_2(area_ids);
        pe_dat_3 = spike_data_table.spike_permEn_3(area_ids);
        pe_dat_4 = spike_data_table.spike_permEn_4(area_ids);
        
        % re-arrange data PE data
        n_neuron_all = length(neuron_ids);
        pe_dat_2 = reshape(pe_dat_2, n_trials_all, n_neuron_all);
        pe_dat_3 = reshape(pe_dat_3, n_trials_all, n_neuron_all);
        pe_dat_4 = reshape(pe_dat_4, n_trials_all, n_neuron_all);

        % drop neurons that are very silent
        pe_dat_2(:, (mean(pe_dat_2))<.0001) = [];
        pe_dat_3(:, (mean(pe_dat_2))<.0001) = [];
        pe_dat_4(:, (mean(pe_dat_2))<.0001) = [];
        
        % also save average PE of remaining neurons
        avg_pe_mat_2 = mean(pe_dat_2);
        avg_pe_mat_3 = mean(pe_dat_3);
        avg_pe_mat_4 = mean(pe_dat_4);

        % another check if we still have neurons left
        if size(pe_dat_2,2) < 2
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
        for i = 1:length(layers)
            % vgg sd 
            col_name = [layers{i} '_layer_nz_sd_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(vgg_pc1_table(:, col_name), ... 
                area_ids, n_trials_all, n_neuron_all, drop_trials);        
            % vgg sum 
            col_name = [layers{i} '_layer_sum_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(vgg_pc1_table(:, col_name), ... 
                area_ids, n_trials_all, n_neuron_all, drop_trials); 
            % vgg non-zero number 
            col_name = [layers{i} '_layer_nz_num_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(vgg_pc1_table(:, col_name), ... 
                area_ids, n_trials_all, n_neuron_all, drop_trials);  
        end
        % prepare hmax data
        layers = {'c1', 'c2'};
        for i = 1:length(layers)
            % hmax sd 
            col_name = [layers{i} '_sd_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(hmax_pc1_table(:, col_name), ... 
                area_ids, n_trials_all, n_neuron_all, drop_trials); 
            % hmax med
            col_name = [layers{i} '_sum_score'];
            stim_data.(col_name(1:end-6)) = prepare_pls_input(hmax_pc1_table(:, col_name), ... 
                area_ids, n_trials_all, n_neuron_all, drop_trials);
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
            stim_data.first_layer_nz_sd,...                                % vgg sd
            stim_data.second_layer_nz_sd,...
            stim_data.third_layer_nz_sd,...
            stim_data.fourth_layer_nz_sd,...
            stim_data.fifth_layer_nz_sd,...
            stim_data.first_layer_nz_num,...                               % vgg num
            stim_data.second_layer_nz_num,...
            stim_data.third_layer_nz_num,...
            stim_data.fourth_layer_nz_num,...
            stim_data.fifth_layer_nz_num,...
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
            stim_data.third_layer_nz_sd,...
            stim_data.fourth_layer_nz_sd,...
            stim_data.fifth_layer_nz_sd,...
            stim_data.third_layer_nz_num,...
            stim_data.fourth_layer_nz_num,...
            stim_data.fifth_layer_nz_num,...
            stim_data.c2_sum,...
            stim_data.c2_sd];
        result_late = pls_analysis({neuro_data}, n_trails_used, 1, option);

        % run one model using early layers
        option.stacked_behavdata = [
            stim_data.first_layer_sum,...                                  % vgg sum
            stim_data.second_layer_sum,...
            stim_data.third_layer_sum,...
            stim_data.first_layer_nz_sd,...                                % vgg sd
            stim_data.second_layer_nz_sd,...
            stim_data.third_layer_nz_sd,...
            stim_data.first_layer_nz_num,...                               % vgg num
            stim_data.second_layer_nz_num,...
            stim_data.third_layer_nz_num,...
            stim_data.c1_sum,...                                           % hmax sum
            stim_data.c1_sd];                                              % hmax sd                                          
        result_early = pls_analysis({neuro_data}, n_trails_used, 1, option);

        % save subject-specific stimulus weights
        behav_weights.all_layers{pat_id,i_ses} = result_combined.v;
        behav_weights.late_layers{pat_id,i_ses} = result_late.v;
        behav_weights.early_layers{pat_id,i_ses} = result_early.v;
        
        neuro_weights.all_layers{pat_id,i_ses} = result_combined.u;
        neuro_weights.late_layers{pat_id,i_ses} = result_late.u;
        neuro_weights.early_layers{pat_id,i_ses} = result_early.u;

        % save number of trials
        trial_nums(pat_id, i_ses) = n_trails_used;
        % number of neurons
        neuron_nums(pat_id, i_ses) = n_neuron_all;

        % save LVLV corrs
        lvlv_corrs.all_layers{pat_id,i_ses} = result_combined.lvlvcorr(1);
        lvlv_corrs.late_layers{pat_id,i_ses} = result_late.lvlvcorr(1);
        lvlv_corrs.early_layers{pat_id,i_ses} = result_early.lvlvcorr(1);

        % also get bounds of LVLV corrs  
        % use the actual lvlv corr distribution
        lvlv_corr_lims.all_layers{pat_id,i_ses} = [percentile(squeeze(...
            result_combined.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_combined.boot_result.lvlvcorrdistrib(1,:)),99)];
        lvlv_corr_lims.late_layers{pat_id,i_ses} = [percentile(squeeze(...
            result_late.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_late.boot_result.lvlvcorrdistrib(1,:)),99)];
        lvlv_corr_lims.early_layers{pat_id,i_ses} = [percentile(squeeze(...
            result_early.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_early.boot_result.lvlvcorrdistrib(1,:)),99)];
    end
    break
end

pls_results = struct(... 
    behav_weights = behav_weights, ... 
    neuro_weights = neuro_weights, ...
    lvlv_corrs = lvlv_corrs, ... 
    lvlv_corr_lims = lvlv_corr_lims, ...
    trial_nums = trial_nums, ...
    neuron_nums = neuron_nums)

% save 
save([pls_dir + 'SpikeVar_PLS_' + region + '_results.mat'], 'pls_results')

end

