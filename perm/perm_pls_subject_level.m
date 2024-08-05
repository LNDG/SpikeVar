function [pls_results, permutations, info] = perm_pls_subject_level(spike_data_table, hmax_pc1_table, vgg_pc1_table)
%pls_subject_level Function that loops over all subjects and performs 3 behavioral pls
%   on for each subject.
%   Takes as input the spike_data_table, and the first pca of hmax and vgg 
%   outputs.

%% initialize
% get IDs
pat_ids = unique(spike_data_table.Participant);
n_pat = length(pat_ids);

n_perm = 1000;
pe_thresh = 0.0001;
trial_thresh = 1/3;

% initialize info
trial_nums = num2cell(zeros(n_pat,1));
drop_trial_nums = num2cell(zeros(n_pat,1));
total_trials = num2cell(zeros(n_pat,1));

neuron_nums = num2cell(zeros(n_pat,1));
drop_neuron_nums = num2cell(zeros(n_pat,1));
total_neurons = num2cell(zeros(n_pat,1)); 

all_neuron_ids = num2cell(zeros(n_pat,1)); 
avg_pe = num2cell(zeros(n_pat,1));

i_ses = 1;

regions = {'hippocampus', 'amygdala'};

%% loop across regions and patients
for region = regions
    display(region{1});   
    count = 1; 
    for pat_id = pat_ids'
        % set of neurons (all, hippocampus <3  or amygdala =>3)
        % 1 = right Hippocampius, 2 = left Hippocampus
        % 3 = right Amygdala, 4 = left Amygdala
        if region{1} == "hippocampus"
            area_ids = find(spike_data_table.Participant==pat_id &...
                           spike_data_table.brainArea < 3 & ... 
                           spike_data_table.Session==i_ses);
        elseif region == "amygdala"
            area_ids = find(spike_data_table.Participant==pat_id &...
                           spike_data_table.brainArea >= 3 & ... 
                           spike_data_table.Session==i_ses);
        else 
            fprintf('Region has to hippocampus or amygdala')
        end
        
        % only continue if we have data
        if isempty(area_ids)
            count = count + 1; 
            continue
        end
     
        % get number of trials & neuron id
        n_trials_all = unique(spike_data_table.NumTrials(spike_data_table.Participant==pat_id &...
             spike_data_table.brainArea<3 & spike_data_table.Session==i_ses));
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
        n_neurons_initial = size(pe_dat_2,2);
        pe_dat_2(:, (mean(pe_dat_2))<=pe_thresh) = [];
        pe_dat_3(:, (mean(pe_dat_3))<=pe_thresh) = [];
        pe_dat_4(:, (mean(pe_dat_4))<=pe_thresh) = [];
        spike_sum_dat(:, (mean(pe_dat_2))<=pe_thresh) = [];
        neuron_ids(mean(pe_dat_2)<=pe_thresh) = [];

        % another check if we still have neurons left
        n_neurons_used = size(pe_dat_2,2);
        n_drop_neurons = n_neurons_initial - n_neurons_used;
 
        % get trials in which less than trial threshold of neurons are
        drop_trials = find(sum(sign(pe_dat_2),2)<=...
            (size(pe_dat_2,2)*trial_thresh));
%         drop_trials_all.(["ses_"+num2str(i_ses)]){pat_id} = drop_trials;
        n_drop_trials = length(drop_trials);
            
        % initialize the spiking input for the pls and drop silent trials
        pe_dat = [pe_dat_2, pe_dat_3, pe_dat_4];
        pe_dat(drop_trials,:) = [];
        spike_sum_dat(drop_trials,:) = [];
        n_trials_used = size(pe_dat,1);
        
        % save average PE of PLS input
        avg_pe{pat_id,1} = mean(pe_dat, "all");
               
        % IMAGE DATA
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
        
        % save number of trials neurons etc
        trial_nums{pat_id,1} = n_trials_used;
        drop_trial_nums{pat_id,1} = n_drop_trials;
        neuron_nums{pat_id,1} = n_neurons_used;
        drop_neuron_nums{pat_id,1} = n_drop_neurons;
        all_neuron_ids{pat_id,1} = neuron_ids;
        total_neurons{pat_id,1} = n_neuron_all; 
        total_trials{pat_id,1} = n_trials_all;

        % run PLS models for PE
        [result_early, result_late, result_all] = run_pls_models(pe_dat, stim_data, n_trials_used);             
        % save LVLV corrs
        pls_results.pe.(region{1}).all_layers(pat_id,1) = result_all.lvlvcorr(1);
        pls_results.pe.(region{1}).late_layers(pat_id,1) = result_late.lvlvcorr(1);
        pls_results.pe.(region{1}).early_layers(pat_id,1) = result_early.lvlvcorr(1);    

        if region{1} == "hippocampus"
            % run PLS models for spike rate
            [result_early, result_late, result_all] = run_pls_models(spike_sum_dat, stim_data, n_trials_used);             
            % save LVLV corrs
            pls_results.spike_rate.hpc.all_layers(pat_id,1) = result_all.lvlvcorr(1);
            pls_results.spike_rate.hpc.late_layers(pat_id,1) = result_late.lvlvcorr(1);
            pls_results.spike_rate.hpc.early_layers(pat_id,1) = result_early.lvlvcorr(1);
        end
                
        % permutations
        for p = 1:n_perm
            % permute dataset
            r = randperm(n_trials_used);
            pe_dat_perm = pe_dat(r,:);
            spike_sum_dat_perm = spike_sum_dat(r,:);
            
            % run PLS models for PE
            [result_early, result_late, result_all] = perm_run_pls_models(pe_dat_perm, stim_data, n_trials_used);             
            % save LVLV corrs
            permutations.pe.(region{1}).all_layers(pat_id,p) = result_all.lvlvcorr(1);
            permutations.pe.(region{1}).late_layers(pat_id,p) = result_late.lvlvcorr(1);
            permutations.pe.(region{1}).early_layers(pat_id,p) = result_early.lvlvcorr(1);    

            if region{1} == "hippocampus"
                % run PLS models for spike rate
                [result_early, result_late, result_all] = perm_run_pls_models(spike_sum_dat_perm, stim_data, n_trials_used);             
                % save LVLV corrs
                permutations.spike_rate.hpc.all_layers(pat_id,p) = result_all.lvlvcorr(1);
                permutations.spike_rate.hpc.late_layers(pat_id,p) = result_late.lvlvcorr(1);
                permutations.spike_rate.hpc.early_layers(pat_id,p) = result_early.lvlvcorr(1);
            end
      
        end % permutation
        count = count + 1; 
    end
    % save info
    info.(region{1}) = struct( ... 
        "pe_thresh",  pe_thresh, ... 
        "trial_thresh", trial_thresh, ...
        "total_trials", total_trials, ... 
        "total_neurons", total_neurons, ... 
        "trial_nums", trial_nums, ... 
        "drop_trial_nums", drop_trial_nums, ... 
        "neuron_nums", neuron_nums, ...
        "drop_neuron_nums", drop_neuron_nums, ...
        "all_neuron_ids", all_neuron_ids, ... 
        "avg_pe", avg_pe);
    % print patient number
    fprintf("\n \n \n Finished PLS for Patient %d, %s. \n", count, region{1})
    pause(.5)
end
end

