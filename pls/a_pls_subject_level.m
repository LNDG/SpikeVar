function pls_results = a_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table, pls_dir, region, measure)
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
        if i_ses == 2
            % for now we only analyse the encoding phase
            continue
        end 
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
                           spike_data_table.Session==i_ses);
        else 
            fprintf('Region has to hippocampus or amygdala')
        end
        
        % only continue if we have data
        if isempty(area_ids)
            continue
        end
     
        % get number of trials (varies between participants and sessions)
        n_trials_all = unique(spike_data_table.NumTrials(spike_data_table.Participant==pat_id &...
             spike_data_table.brainArea<3 & spike_data_table.Session==i_ses));
        % save neuron info
        neuron_ids = unique(spike_data_table.NeuronCount(area_ids));
        n_neuron_all = length(neuron_ids);
        
        % extract spike sum and re-arrange data
        spike_sum_dat = sum(spike_data_table.stimuliSpikeCounts(area_ids,1:100),2);
        spike_sum_dat = reshape(spike_sum_dat , n_trials_all, n_neuron_all);

        % extract PE for different motif lengths
        pe_dat_2 = spike_data_table.spike_permEn_2(area_ids);
        pe_dat_3 = spike_data_table.spike_permEn_3(area_ids);
        pe_dat_4 = spike_data_table.spike_permEn_4(area_ids);
        
        % re-arrange data PE data
        pe_dat_2 = reshape(pe_dat_2, n_trials_all, n_neuron_all);
        pe_dat_3 = reshape(pe_dat_3, n_trials_all, n_neuron_all);
        pe_dat_4 = reshape(pe_dat_4, n_trials_all, n_neuron_all);

        % drop neurons that are very silent
        pe_dat_2(:, (mean(pe_dat_2))<.0001) = [];
        pe_dat_3(:, (mean(pe_dat_3))<.0001) = [];
        pe_dat_4(:, (mean(pe_dat_4))<.0001) = [];
        spike_sum_dat(:, (mean(pe_dat_2))<.0001) = [];
        neuron_ids(mean(pe_dat_2)<.0001) = [];
        
        % also save average PE of remaining neurons
        avg_pe_mat_2 = mean(pe_dat_2);
        avg_pe_mat_3 = mean(pe_dat_3);
        avg_pe_mat_4 = mean(pe_dat_4);

        % another check if we still have neurons left
        n_neurons_used = size(pe_dat_2,2);
        %if size(pe_dat_2,2) < 2
        %    continue
        %end
               
        % get trials in which less than 1/3 of neurons are
        % spiking at all
        drop_trials = find(sum(sign(pe_dat_2),2)<=...
            (size(pe_dat_2,2)/3));        
        
        % initialize the spiking input for the pls and drop silent trials
        pe_dat = [pe_dat_2, pe_dat_3, pe_dat_4];
        pe_dat(drop_trials,:) = [];
        spike_sum_dat(drop_trials,:) = [];
        n_trials_used = size(pe_dat,1);
        
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
        if measure == "SpikeNum"
            [result_early, result_late, result_all] = run_pls_models(spike_sum_dat, stim_data, n_trials_used);
        elseif measure == "PE"
            [result_early, result_late, result_all] = run_pls_models(pe_dat, stim_data, n_trials_used);
        end
        
        % save subject-specific stimulus weights
        stim_weights.all_layers{pat_id,i_ses} = result_all.v;
        stim_weights.late_layers{pat_id,i_ses} = result_late.v;
        stim_weights.early_layers{pat_id,i_ses} = result_early.v;
        
        neuro_weights.all_layers{pat_id,i_ses} = result_all.u;
        neuro_weights.late_layers{pat_id,i_ses} = result_late.u;
        neuro_weights.early_layers{pat_id,i_ses} = result_early.u;

        % save number of trials
        trial_nums(pat_id, i_ses) = n_trials_used;
        % number of neurons
        neuron_nums(pat_id, i_ses) = n_neurons_used;
        % all neuron ids
        all_neuron_ids{pat_id, i_ses} = neuron_ids;

        % save LVLV corrs
        lvlv_corrs.all_layers{i_ses}(pat_id) = result_all.lvlvcorr(1);
        lvlv_corrs.late_layers{i_ses}(pat_id) = result_late.lvlvcorr(1);
        lvlv_corrs.early_layers{i_ses}(pat_id) = result_early.lvlvcorr(1);

        % also get bounds of LVLV corrs  
        % use the actual lvlv corr distribution
        lvlv_corr_lims.all_layers{i_ses}(pat_id,:) = [percentile(squeeze(...
            result_all.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_all.boot_result.lvlvcorrdistrib(1,:)),99)];
        lvlv_corr_lims.late_layers{i_ses}(pat_id,:) = [percentile(squeeze(...
            result_late.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_late.boot_result.lvlvcorrdistrib(1,:)),99)];
        lvlv_corr_lims.early_layers{i_ses}(pat_id,:) = [percentile(squeeze(...
            result_early.boot_result.lvlvcorrdistrib(1,:)),1),...
            percentile(squeeze(...
            result_early.boot_result.lvlvcorrdistrib(1,:)),99)];
    end
    fprintf("\n \n \n Finished PLS for Patient %d.", count)
    count = count + 1;
end

pls_results = struct(... 
    stim_weights = stim_weights, ... 
    neuro_weights = neuro_weights, ...
    lvlv_corrs = lvlv_corrs, ... 
    lvlv_corr_lims = lvlv_corr_lims); 

pls_results.trial_nums = trial_nums; 
pls_results.neuron_nums = neuron_nums;
pls_results.all_neuron_ids = all_neuron_ids;

% save 
save([pls_dir + 'SpikeVar_' + measure + '_PLS_' + region + '_results.mat'], ... 
    'pls_results')

end

