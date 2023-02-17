function save_abs_stim_weights(stim_weights, relevant_ids, pls_dir)
% Get the absoulute stimulus weight for each layer
    for sub_id = 1:size(stim_weights.late_layers,1)
        if ~isempty(stim_weights.late_layers{sub_id,1})
            all_abs_stim_weights_late(:,sub_id) = (abs(stim_weights.late_layers{sub_id,1}(:,1)));
        else
            all_abs_stim_weights_late(:,sub_id) = NaN(size(stim_weights.late_layers{1,1},1),1);
        end
    end

    for lay = 1:size(stim_weights.late_layers{1},1)
        abs_stim_weights_late{lay,1} = (all_abs_stim_weights_late(lay,relevant_ids))';
    end

    % save weights to run repeated measure models else where
    vggweight_rm_tab = array2table((all_abs_stim_weights_late(3:8,relevant_ids))');
    vggweight_rm_tab.Properties.VariableNames = {'VGGsumweight1', 'VGGsumweight2',...
        'VGGsumweight3', 'VGGSDweight1', 'VGGSDweight2', 'VGGSDweight3'};
    writetable(vggweight_rm_tab, [pls_dir + 'SpikeVar_VGGweight_table.csv'])
end