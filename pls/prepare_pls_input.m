function pca_table = prepare_pls_input(pca_table, area_ids, num_trials, n_neurons, drop_trials)
% prepare the input for the pls: Take only ids of relevante brain area;
% reshape to matrix; take only one column; drop bad trials
    pca_table = pca_table(area_ids, :);
    pca_table = table2array(pca_table);
    pca_table = reshape(pca_table, num_trials,...
        n_neurons);
    pca_table = squeeze(pca_table(:,1));
    pca_table(drop_trials,:) = [];  
end
