function results = get_layer_wise_pca(layer_wise_estim)
    [results.coeff, results.score, results.latent,~,results.exp_var,~] = ...
        pca(zscore(layer_wise_estim), 'Centered', true,'VariableWeights', 'Variance');

    results.loadings = (corr(results.score,...
        layer_wise_estim(:,1:size(results.score,2))))';
    results.exp_var_1st_pc = results.exp_var(1);
    results.exp_normvar_1st_pc = results.exp_var(1)/...
                        (100/size(results.exp_var,1));
end

