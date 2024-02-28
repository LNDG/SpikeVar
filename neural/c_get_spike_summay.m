function spike_summary = c_get_spike_summay(spike_data_table, neuro_dir)
    pe_table = spike_data_table;
    pat_ids = unique(pe_table.Participant);
    % get task variant codes
    for i_pat = 1:length(pat_ids)
        for ises = 1:2
            if ~isempty(find(pe_table.Participant==...
                    pat_ids(i_pat)& pe_table.Session== ises))
                rel_it = find(pe_table.Participant==...
                    pat_ids(i_pat)& pe_table.Session== ises);
                task_vars(i_pat,ises) = unique(pe_table.TaskVariant(rel_it));
                image_dur(i_pat,ises) = unique(pe_table.learnImageDur(rel_it))/1000 ;
                avg_pe(i_pat,ises) = mean(pe_table.spike_permEn_4(rel_it));
                spike_rate(i_pat,ises) =  mean(sum(pe_table.stimuliSpikeCounts(rel_it,1:image_dur(i_pat,ises)*100-1),2))/image_dur(i_pat,ises);
            else
                task_vars(i_pat,ises) = NaN;
                image_dur(i_pat,ises) = NaN;
                avg_pe(i_pat,ises) = NaN;
                spike_rate(i_pat,ises) = NaN;
            end
        end
    end
    % binarize task variant
    task_vars(task_vars<3) = 0;
    task_vars(task_vars ==3) = 1;
    spike_summary = struct('task_vars', task_vars, ... 
        'image_dur', image_dur, ...
        'avg_pe', avg_pe, ...
        'spike_rate', spike_rate);
    save([neuro_dir + 'SpikeVar_spike_summary.mat'], 'spike_summary')
end