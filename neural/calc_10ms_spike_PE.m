clear
clc
% check if we are in the correct directory, change it if needed
base_dir = '/Users/kamp/PhD/spikevar';
neuro_dir = fullfile(base_dir, 'repo', 'SpikeVar', 'neural/');
data_dir = fullfile(base_dir,'data', 'neuro/');
stim_dir = fullfile(base_dir,'stimuli_rn/');
tool_dir = fullfile(base_dir, 'repo', 'SpikeVar', 'toolboxes');

addpath(tool_dir);
load([data_dir 'SpikeVar_spike_data.mat'], 'spike_data_table');
%%
% get subject indices 
pat_ids = unique(spike_data_table.Participant);


for participant_id = pat_ids'
    % get session indices
    ses_ids = unique(spike_data_table.Session(spike_data_table.Participant==participant_id));
    for session_id = ses_ids'
        % get neuron indices
        neuron_ids = unique(spike_data_table.Neuron(spike_data_table.Participant==participant_id &...
            spike_data_table.Session==session_id));
        for neuron_id = neuron_ids'
            % get trials indices
            trial_ids = find(spike_data_table.Participant==participant_id & ...
                spike_data_table.Session==session_id & spike_data_table.Neuron==neuron_id);
            % loop over trials and extract PE
            for trial_id = trial_ids'
                % only use the first 1000 ms of data after stim onset
                % in the case of 10 ms bins this results in 100 bins
                stimuliSpikeCounts = [spike_data_table.stimuliSpikeCounts(trial_id,1:100)];

                % save number of spikes per trial and neuron
                spike_data_table.numSpikes(trial_id) = ...
                    nansum(stimuliSpikeCounts);

                % get permEn for m = 2-4
                [tt,perm]=PermEn(stimuliSpikeCounts, 'm', 4, 'tau',1,...
                    'Norm',true);
                spike_data_table.spike_permEn_2(trial_id) = ...
                    (perm(2));
                spike_data_table.spike_permEn_3(trial_id) = ...
                    (perm(3));
                spike_data_table.spike_permEn_4(trial_id) = ...
                    (perm(4));
            end
        end
    end
    disp(['Done with participant ' num2str(participant_id)])
end
%% save table inlcuding single trial PE
save([data_dir 'SpikeVar_spike_data_PE.mat'], 'spike_data_table')
