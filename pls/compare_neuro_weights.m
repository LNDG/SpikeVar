function compare_neuro_weights(pls_results, VS_cells, MS_cells, rel_ids, all_neuron_ids)
    %% compare neural weights  between neurons
    % extract info per subject
    for sub_i = 1:length(rel_ids)
        % extract neuron infos
       n_neurons = length(all_neuron_ids{rel_ids(sub_i),1})
       ms_neu_inds{rel_ids(sub_i), 1} = find(ismember(all_neuron_ids{rel_ids(sub_i),1},...
           [MS_cells MS_cells+n_neurons(rel_ids(sub_i), 1),...
           MS_cells+(n_neurons(rel_ids(sub_i), 1))*2]));
       vs_neu_inds{rel_ids(sub_i), 1} = find(ismember(all_neuron_ids{rel_ids(sub_i),1},...
           [VS_cells VS_cells+n_neurons(rel_ids(sub_i), 1),...
           VS_cells+(n_neurons(rel_ids(sub_i), 1))*2]));
       % extract weights
       ms_neu_wghts{rel_ids(sub_i), 1} = abs(pls_results.neuro_weights.late_layers{rel_ids(sub_i), 1}(...
           ms_neu_inds{rel_ids(sub_i), 1}));
       vs_neu_wghts{rel_ids(sub_i), 1} = abs(pls_results.neuro_weights.late_layers{rel_ids(sub_i), 1}(...
           vs_neu_inds{rel_ids(sub_i), 1}));
    end

    % plot comparison
    all_ms_neu_wghts = vertcat(ms_neu_wghts{:});
    all_vs_neu_wghts = vertcat(vs_neu_wghts{:});

    figure
    histogram(all_vs_neu_wghts)
    hold on
    histogram(all_ms_neu_wghts)
    legend('Vision selective neurons', 'Memory selective neurons')
    xlabel('abs(weight)')
    ylabel('# occurrence')

    % compare both distributions, use non-parametric
    [p,h,stats] = ranksum(all_vs_neu_wghts, all_ms_neu_wghts);
end