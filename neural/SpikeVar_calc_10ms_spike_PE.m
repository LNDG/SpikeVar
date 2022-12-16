clear
clc
% check if we are in the correct directory, change it if needed
cd('/Users/waschke/OneDrive/Old data and code/SV/OurCode/')
% load table with all relevant info,
load('SPSSinput_10_msec.mat')
% add essential scripts
addpath('/Users/waschke/Documents/MATLAB/Code/Scripts/Essentials/')
addpath(genpath('/Users/waschke/Documents/Matlabtoolboxes/EntropyHub'))
%%
% check if table inlcuding PermEN already exists, no need to redo in that
% case
checko = 1;%exist('SpikeVar_all_trials_incl_10_ms_PE.mat');

if checko ~=2
    disp('Calculating PE per neuron and trial. This might take some time...')
    % use all neurons
    perm_active_table = fulltable;
    cur_pat_ids = unique(perm_active_table.Participant);
    
    for pid = 1:length(cur_pat_ids)
        cid = cur_pat_ids(pid);
        % got to bin on a neuron by neuron basis
        % neuron numbers
        for ises = 1:length(unique(perm_active_table.Session(perm_active_table.Participant==cid)))
            
            neu_nums =unique(perm_active_table.Neuron(perm_active_table.Participant==cid &...
                perm_active_table.Session==ises));
            for in = 1:length(neu_nums)
                
                % determine number of trials within current participant
                numts = find(perm_active_table.Participant==cid &...
                    perm_active_table.Session==ises & perm_active_table.Neuron==neu_nums(in));
                % loop over trials and extract PE
                for its=1:length(numts)

                    % only use the first 1000 ms of data after stim onset
                    % In the case of 10 ms bins this results in 100 bins
                    rel_psth10 = [perm_active_table.psth10(numts(its),1:100)];
                    
                    % save number of spikes per trial and neuron
                    perm_active_table.numspikes10(numts(its)) = ...
                        nansum(rel_psth10);
                    
                    % get permEn for m = 2-4
                    [tt,temp_perm]=PermEn(rel_psth10, 'm', 4, 'tau',1,...
                        'Norm',true);
                    perm_active_table.psth10_permEN_2(numts(its)) = ...
                        (temp_perm(2));
                    perm_active_table.psth10_permEN_3(numts(its)) = ...
                        (temp_perm(3));
                    perm_active_table.psth10_permEN_4(numts(its)) = ...
                        (temp_perm(4));
                end
            end
        end
        disp(['Done with participant ' num2str(pid)])
    end
    %% save table inlcuding single trial PE
    save('SpikeVar_all_trials_incl_10_ms_PE.mat', 'perm_active_table', '-v7.3')
else
    disp('PE already calculated, loading file...')
    load('SpikeVar_all_trials_incl_10_ms_PE.mat')
    cur_pat_ids = unique(perm_active_table.Participant);
end