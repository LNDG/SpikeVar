function results = apply_pca(feature_wise_estim)
    
    % check for features without variance
    feature_var = var(feature_wise_estim, 0, 1); 
    
    if ~ all(feature_var)  
        feature_wise_estim(:, feature_var==0) = [];
    fi 
    
    [results.coeff,results.score,...
     results.latent,~,results.exp_var,~] = ...
        pca(zscore(feature_wise_estim), 'Centered', true,'VariableWeights', 'Variance');

    results.loadings = (corr(results.score,...
        feature_wise_estim(:,1:size(results.score,2))))';
    results.exp_var_1st_pc = results.exp_var(1);
    results.exp_normvar_1st_pc = results.exp_var(1)/...
                        (100/size(results.exp_var,1));
end

